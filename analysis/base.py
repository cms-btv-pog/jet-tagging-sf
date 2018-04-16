# -*- coding: utf-8 -*-


__all__ = ["AnalysisTask", "DatasetTask", "GridWorkflow"]


import re
import os

import law
import luigi

from analysis.config_2017 import analysis as analysis_2017, campaign as campaign_2017


law.contrib.load("cms", "git", "glite", "tasks", "wlcg")


class AnalysisTask(law.Task):

    version = luigi.Parameter()

    outputs_siblings = True

    def __init__(self, *args, **kwargs):
        super(AnalysisTask, self).__init__(*args, **kwargs)

        self.analysis_inst = analysis_2017
        self.campaign_inst = campaign_2017

        self.config_inst = self.analysis_inst.get_config(self.campaign_inst.id)

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


class UploadCMSSW(AnalysisTask, law.BundleCMSSW, law.TransferLocalFile):

    # settings for BunddleCMSSW
    cmssw_path = os.getenv("CMSSW_BASE")

    # settings for TransferLocalFile
    source_path = None

    version = None
    task_namespace = None

    def single_output(self):
        path = "{}.tgz".format(os.path.basename(self.cmssw_path))
        return self.wlcg_target(path)

    def output(self):
        return law.TransferLocalFile.output(self)

    def run(self):
        bundle = law.LocalFileTarget(is_tmp="tgz")
        self.bundle(bundle)
        self.transfer(bundle)


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
