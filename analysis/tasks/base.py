# -*- coding: utf-8 -*-


__all__ = ["AnalysisTask", "DatasetTask", "GridWorkflow"]


import re
import os
import shutil

import law
import luigi

from analysis.config.jet_tagging_sf import analysis
from analysis.util import calc_checksum


law.contrib.load("cms", "git", "glite", "tasks", "wlcg")


class AnalysisTask(law.Task):

    version = luigi.Parameter()

    outputs_siblings = True

    # the campaign is hardcoded for the moment
    campaign = "2017_Run2_pp_13TeV_ICHEP18"

    def __init__(self, *args, **kwargs):
        super(AnalysisTask, self).__init__(*args, **kwargs)

        self.analysis_inst = analysis
        self.config_inst = self.analysis_inst.get_config(self.campaign)

    def store_parts(self):
        parts = (self.__class__.__name__, self.config_inst.name)
        if self.version is not None:
            parts += (self.version,)
        return parts

    def local_path(self, *path):
        return os.path.join(os.environ["JTSF_STORE"], *(self.store_parts() + path))

    def wlcg_path(self, *path):
        return os.path.join(*(self.store_parts() + path))

    def local_target(self, *args):
        cls = law.LocalFileTarget if args else law.LocalDirectoryTarget
        return cls(self.local_path(*args))

    def wlcg_target(self, *args, **kwargs):
        cls = law.WLCGFileTarget if args else law.WLCGDirectoryTarget
        return cls(self.wlcg_path(*args), **kwargs)


class DatasetTask(AnalysisTask):

    dataset = luigi.Parameter(default="data_ee")

    def __init__(self, *args, **kwargs):
        super(DatasetTask, self).__init__(*args, **kwargs)

        self.dataset_inst = self.config_inst.get_dataset(self.dataset)

    def store_parts(self):
        return super(DatasetTask, self).store_parts() + (self.dataset,)

    def create_branch_map(self):
        return list(range(self.dataset_inst.n_files))

    def glite_output_postfix(self):
        self.get_branch_map()
        return "_{}To{}".format(self.start_branch, self.end_branch)


class GridWorkflow(AnalysisTask, law.GLiteWorkflow):

    glite_ce_map = {
        "RWTH": "grid-ce.physik.rwth-aachen.de:8443/cream-pbs-cms",
        "RWTH_short": "grid-ce.physik.rwth-aachen.de:8443/cream-pbs-short",
        "CNAF": [
            "ce04-lcg.cr.cnaf.infn.it:8443/cream-lsf-cms",
            "ce05-lcg.cr.cnaf.infn.it:8443/cream-lsf-cms",
            "ce06-lcg.cr.cnaf.infn.it:8443/cream-lsf-cms",
        ],
    }

    grid_ce = law.CSVParameter(default=["RWTH"], significant=False, description="target computing "
        "element(s)")

    @classmethod
    def modify_param_values(cls, params):
        params = AnalysisTask.modify_param_values(params)
        if "workflow" in params and law.is_no_param(params["workflow"]):
            grid_ce = params["grid_ce"]
            ces = []
            for ce in grid_ce:
                ces.append(cls.glite_ce_map.get(ce, ce))
            params["glite_ce"] = tuple(law.util.flatten(ces))
        return params

    def glite_workflow_requires(self):
        reqs = law.GLiteWorkflow.glite_workflow_requires(self)
        reqs["cmssw"] = UploadCMSSW.req(self, replicas=10, _prefer_cli=["replicas"])
        reqs["software"] = UploadSoftware.req(self, replicas=10, _prefer_cli=["replicas"])
        reqs["repo"] = UploadRepo.req(self, replicas=10, _prefer_cli=["replicas"])
        return reqs

    def glite_output_directory(self):
        return law.WLCGDirectoryTarget(self.wlcg_path())

    def glite_output_uri(self):
        return self.glite_output_directory().url(cmd="listdir")

    def glite_bootstrap_file(self):
        return law.util.rel_path(__file__, "grid_bootstrap.sh")

    def glite_job_config(self, config, job_num, branches):
        config = law.GLiteWorkflow.glite_job_config(self, config, job_num, branches)
        reqs = self.glite_workflow_requires()
        config.vo = "cms:/cms/dcms"
        config.render_variables["cmssw_base"] = reqs["cmssw"].output().dir.url()
        config.render_variables["cmssw_version"] = os.getenv("CMSSW_VERSION")
        config.render_variables["software_base"] = reqs["software"].output().dir.url()
        config.render_variables["repo_checksum"] = reqs["repo"].checksum
        config.render_variables["repo_base"] = reqs["repo"].output().dir.url()
        return config


class InstallCMSSWCode(AnalysisTask):

    version = None

    def __init__(self, *args, **kwargs):
        super(InstallCMSSWCode, self).__init__(*args, **kwargs)

        self._checksum = None

    @property
    def checksum(self):
        if self._checksum is None:
            path = os.path.join(os.getenv("JTSF_BASE"), "cmssw")
            self._checksum = calc_checksum(path, exclude=["*.pyc", "*.git*", "tmpfiles*"])

        return self._checksum

    def output(self):
        return self.local_target("{}.txt".format(self.checksum))

    def run(self):
        # copy the current cmssw code to the CMSSW_BASE directory
        for subsystem in ["jet_tagging_sf"]:
            src = os.path.join(os.getenv("JTSF_BASE"), "cmssw", subsystem)
            dst = os.path.join(os.getenv("CMSSW_BASE"), "src", subsystem)
            if os.path.exists(dst):
                shutil.rmtree(dst)
            shutil.copytree(src, dst)

        # install the software
        code = law.util.interruptable_popen("scram b", shell=True, executable="/bin/bash",
            cwd=os.path.dirname(dst))[0]
        if code != 0:
            raise Exception("scram build failed")

        # touch the flag output file
        output = self.output()
        output.parent.touch(0o0770)
        output.touch(self.checksum)


class UploadCMSSW(AnalysisTask, law.BundleCMSSW, law.TransferLocalFile):

    force_upload = luigi.BoolParameter(default=False, description="force uploading")

    # settings for BunddleCMSSW
    cmssw_path = os.getenv("CMSSW_BASE")

    # settings for TransferLocalFile
    source_path = None

    version = None
    task_namespace = None

    def __init__(self, *args, **kwargs):
        super(UploadCMSSW, self).__init__(*args, **kwargs)

        self.has_run = False

    def complete(self):
        if self.force_upload:
            return self.has_run
        else:
            return super(UploadCMSSW, self).complete()

    def single_output(self):
        path = "{}.tgz".format(os.path.basename(self.cmssw_path))
        return self.wlcg_target(path)

    def output(self):
        return law.TransferLocalFile.output(self)

    def run(self):
        bundle = law.LocalFileTarget(is_tmp="tgz")
        self.bundle(bundle)
        self.transfer(bundle)

        self.has_run = True


class UploadSoftware(AnalysisTask, law.TransferLocalFile):

    version = None

    source_path = os.environ["JTSF_SOFTWARE"] + ".tgz"

    def single_output(self):
        return self.wlcg_target(os.path.basename(self.source_path))

    def run(self):
        # create the local bundle
        bundle = law.LocalFileTarget(self.source_path, is_tmp=True)
        def _filter(tarinfo):
            return None if re.search("(\.pyc|\/\.git|\.tgz)$", tarinfo.name) else tarinfo
        bundle.dump(os.path.splitext(self.source_path)[0], filter=_filter)
        self.publish_message("bundled software archive")

        # super run will upload all files for us
        super(UploadSoftware, self).run()


class UploadRepo(AnalysisTask, law.BundleGitRepository, law.TransferLocalFile):

    # settings for BundleGitRepository
    repo_path = os.environ["JTSF_BASE"]

    # settings for TransferLocalFile
    source_path = None

    version = None
    task_namespace = None

    def single_output(self):
        path = "{}.{}.tgz".format(os.path.basename(self.repo_path), self.checksum)
        return self.wlcg_target(path)

    def output(self):
        return law.TransferLocalFile.output(self)

    def run(self):
        bundle = law.LocalFileTarget(is_tmp="tgz")
        self.bundle(bundle)
        self.transfer(bundle)
