# -*- coding: utf-8 -*-


__all__ = ["AnalysisTask", "DatasetTask", "ShiftTask", "GridWorkflow"]


import re
import os
import abc
import shutil
import random
import itertools
import collections
import subprocess

import law
import luigi
import six

from law.workflow.base import BaseWorkflow
from law.util import make_list, tmp_file, interruptable_popen

from analysis.config.jet_tagging_sf import analysis
from analysis.util import calc_checksum


law.contrib.load("arc", "cms", "git", "glite", "numpy", "tasks", "root", "slack", "wlcg", "htcondor",
    "singularity")


class AnalysisTask(law.Task):

    version = luigi.Parameter()
    notify = law.slack.NotifySlackParameter()

    output_collection_cls = law.target.collection.SiblingFileCollection

    config = luigi.Parameter(default=None)

    accepts_messages = True
    message_cache_size = 20

    exclude_params_req = {"notify"}

    def __init__(self, *args, **kwargs):
        super(AnalysisTask, self).__init__(*args, **kwargs)

        if self.config is None:
            self.config = os.environ["JTSF_CAMPAIGN"]

        self.analysis_inst = analysis
        self.config_inst = self.analysis_inst.get_config(self.config)

    def get_version(self, task_cls):
        family = task_cls if isinstance(task_cls, six.string_types) else task_cls.get_task_family()
        return self.config_inst.get_aux("versions")[family]

    def store_parts(self):
        parts = (self.__class__.__name__, self.config_inst.name)
        if self.version is not None:
            parts += (self.version,)
        return parts

    @classmethod
    def modify_param_values(cls, params):
        return params

    def local_path(self, *path):
        parts = [str(p) for p in self.store_parts() + path]
        return os.path.join(os.environ["JTSF_STORE"], *parts)

    def wlcg_path(self, *path):
        parts = [str(p) for p in self.store_parts() + path]
        return os.path.join(*parts)

    def local_target(self, *args):
        cls = law.LocalFileTarget if args else law.LocalDirectoryTarget
        return cls(self.local_path(*args))

    def wlcg_target(self, *args, **kwargs):
        cls = law.wlcg.WLCGFileTarget if args else law.wlcg.WLCGDirectoryTarget
        return cls(self.wlcg_path(*args), **kwargs)


class DatasetTask(AnalysisTask):

    dataset = luigi.Parameter(default="data_B_ee")

    file_merging = None

    def __init__(self, *args, **kwargs):
        super(DatasetTask, self).__init__(*args, **kwargs)

        self.dataset_inst = self.config_inst.get_dataset(self.dataset)

    def store_parts(self):
        return super(DatasetTask, self).store_parts() + (self.dataset,)

    def create_branch_map(self):
        merging_info = self.config_inst.get_aux("file_merging")
        if isinstance(self.file_merging, six.integer_types):
            n = self.file_merging
        elif self.file_merging in merging_info:
            n = merging_info[self.file_merging].get(self.dataset_inst.name, 1)
        else:
            n = self.dataset_inst.n_files
        return list(range(n))

    def glite_output_postfix(self):
        self.get_branch_map()
        return "_{}To{}".format(self.start_branch, self.end_branch)


class ShiftTask(AnalysisTask):

    shift = luigi.Parameter(default="nominal", significant=False,
        description="systematic shift to apply, default: nominal")
    effective_shift = luigi.Parameter(default=None)

    shifts = set()

    exclude_params_req = {"effective_shift"}

    def __init__(self, *args, **kwargs):
        super(ShiftTask, self).__init__(*args, **kwargs)

    def store_parts(self):
        return super(ShiftTask, self).store_parts() + (self.effective_shift,)

    @classmethod
    def modify_param_values(cls, params):
        params = super(ShiftTask, cls).modify_param_values(params)
        params = cls.get_effective_shift(params)
        return params

    @classmethod
    def get_effective_shift(cls, params):
        if "shift" not in params:
            return params
        # if the task does not implement the provided shift, use the nominal one
        if params["shift"] in cls.shifts:
            params["effective_shift"] = params["shift"]
        else:
            params["effective_shift"] = "nominal"
        return params


class WrapperTask(AnalysisTask, law.WrapperTask):

    datasets = law.CSVParameter(default=[], description="datasets to require")
    shifts = law.CSVParameter(default=[], description="shifts to require")
    skip_datasets = law.CSVParameter(default=[], description="datasets to skip, supports patterns")
    skip_shifts = law.CSVParameter(default=[], description="shifts to skip, supports patterns")
    grid_ces = law.CSVParameter(default=[], description="grid CEs to submit to, chosen randomly")

    exclude_db = True

    def __init__(self, *args, **kwargs):
        super(WrapperTask, self).__init__(*args, **kwargs)

        if not self.datasets:
            self.datasets = self.get_default_datasets()

        if not self.shifts:
            self.shifts = self.get_default_shifts()

        if self.skip_datasets:
            filter_fn = lambda d: not law.util.multi_match(d, self.skip_datasets)
            self.datasets = filter(filter_fn, self.datasets)
        if self.skip_shifts:
            filter_fn = lambda d: not law.util.multi_match(d, self.skip_shifts)
            self.shifts = filter(filter_fn, self.shifts)

    @abc.abstractproperty
    def wrapped_task(self):
        return

    def get_default_datasets(self):
        if issubclass(self.wrapped_task, DatasetTask):
            return [dataset.name for dataset in self.config_inst.datasets]
        else:
            return [None]

    def get_default_shifts(self):
        if issubclass(self.wrapped_task, ShiftTask):
            return self.wrapped_task.shifts
        else:
            return [None]

    def requires(self):
        cls = self.wrapped_task

        def req(dataset, shift):
            kwargs = {}
            if dataset is not None:
                kwargs["dataset"] = dataset
            if shift is not None:
                kwargs["shift"] = shift

            if issubclass(cls, GridWorkflow) and self.grid_ces:
                kwargs["grid_ce"] = [random.choice(self.grid_ces)]
                kwargs["_prefer_cli"] = ["grid_ce"]

            return cls.req(self, **kwargs)

        # get parameters, require shifts only for MC
        params_list = []
        for dataset in self.datasets:
            for shift in self.shifts:
                if dataset is not None and shift is not None:
                    # require shifts only for MC
                    if self.config_inst.get_dataset(dataset).is_data and shift != "nominal":
                        continue
                params_list.append((dataset, shift))

        return collections.OrderedDict([(params, req(*params)) for params in params_list])


class GridWorkflow(AnalysisTask, law.glite.GLiteWorkflow, law.arc.ARCWorkflow):

    glite_ce_map = {
        "RWTH": "grid-ce.physik.rwth-aachen.de:8443/cream-pbs-cms",
        "RWTH_short": "grid-ce.physik.rwth-aachen.de:8443/cream-pbs-short",
        "CNAF": [
            "ce04-lcg.cr.cnaf.infn.it:8443/cream-lsf-cms",
            "ce05-lcg.cr.cnaf.infn.it:8443/cream-lsf-cms",
            "ce06-lcg.cr.cnaf.infn.it:8443/cream-lsf-cms",
        ],
        "IRFU": "node74.datagrid.cea.fr:8443/cream-pbs-cms",
        "IIHE": "cream02.iihe.ac.be:8443/cream-pbs-cms",
        "CIEMAT": [
            "creamce02.ciemat.es:8443/cream-pbs-medium",
            "creamce03.ciemat.es:8443/cream-pbs-medium",
        ],
    }
    arc_ce_map = {
        "DESY": "grid-arcce0.desy.de",
        "KIT": ["arc-{}-kit.gridka.de".format(i) for i in range(1, 6 + 1)],
    }

    sl_distribution_map = collections.defaultdict(lambda: "slc7", {"RWTH": "slc6"})
    req_sandbox = "NO_SANDBOX"

    grid_ce = law.CSVParameter(default=["RWTH"], significant=False, description="target computing "
        "element(s)")

    exclude_params_branch = {"grid_ce"}

    @classmethod
    def modify_param_values(cls, params):
        params = AnalysisTask.modify_param_values(params)
        if "workflow" in params and law.is_no_param(params["workflow"]):
            grid_ce = params["grid_ce"]
            workflow = "arc" if grid_ce[0] in cls.arc_ce_map else "glite"
            ces = []
            for ce in grid_ce:
                ces.append(getattr(cls, workflow + "_ce_map").get(ce, ce))
            params[workflow + "_ce"] = tuple(law.util.flatten(ces))
            params["workflow"] = workflow
        return params

    def _setup_workflow_requires(self, reqs):
        # figure out if the CE runs the same operating system as we are
        # if not, upload the software and cmssw from an appropriate sandbox
        self.sl_dist_version = os.getenv("JTSF_DIST_VERSION")

        if not len(set([self.sl_distribution_map[ce] for ce in self.grid_ce])) == 1:
            raise Exception("Cannot submit to multiple CEs running different distributions.")

        if self.sl_distribution_map[self.grid_ce[0]] != self.sl_dist_version:
            self.req_sandbox = self.sl_distribution_map[self.grid_ce[0]]

        reqs["cmssw"] = UploadCMSSW.req(self, replicas=10, _prefer_cli=["replicas"],
            sandbox=self.config_inst.get_aux("sandboxes")[self.req_sandbox])
        reqs["software"] = UploadSoftware.req(self, replicas=10, _prefer_cli=["replicas"],
            sandbox=self.config_inst.get_aux("sandboxes")[self.req_sandbox])
        reqs["repo"] = UploadRepo.req(self, replicas=10, _prefer_cli=["replicas"])

    def _setup_render_variables(self, config, reqs):
        config.render_variables["jtsf_grid_user"] = os.getenv("JTSF_GRID_USER")
        config.render_variables["jtsf_cmssw_setup"] = os.getenv("JTSF_CMSSW_SETUP")
        config.render_variables["cmssw_base_url"] = reqs["cmssw"].output().dir.uri()

        scram_arch = os.getenv("SCRAM_ARCH")
        if self.req_sandbox != "NO_SANDBOX":
            scram_arch = scram_arch.replace(self.sl_dist_version, self.req_sandbox)
        config.render_variables["scram_arch"] = scram_arch

        config.render_variables["cmssw_version"] = os.getenv("CMSSW_VERSION")
        config.render_variables["software_base_url"] = reqs["software"].output().dir.uri()
        config.render_variables["repo_checksum"] = reqs["repo"].checksum
        config.render_variables["repo_base"] = reqs["repo"].output().dir.uri()

    def glite_workflow_requires(self):
        reqs = law.glite.GLiteWorkflow.glite_workflow_requires(self)
        self._setup_workflow_requires(reqs)
        return reqs

    def glite_output_directory(self):
        return law.wlcg.WLCGDirectoryTarget(self.wlcg_path())

    def glite_output_uri(self):
        return self.glite_output_directory().uri(cmd="listdir")

    def glite_bootstrap_file(self):
        return law.util.rel_path(__file__, "files", "grid_bootstrap.sh")

    def glite_job_config(self, config, job_num, branches):
        config = law.glite.GLiteWorkflow.glite_job_config(self, config, job_num, branches)
        self._setup_render_variables(config, self.glite_workflow_requires())
        config.vo = "cms:/cms/dcms"
        return config

    def arc_workflow_requires(self):
        reqs = law.ARCWorkflow.arc_workflow_requires(self)
        self._setup_workflow_requires(reqs)
        return reqs

    def arc_output_directory(self):
        return self.glite_output_directory()

    def arc_output_uri(self):
        return self.glite_output_uri()

    def arc_bootstrap_file(self):
        return self.glite_bootstrap_file()

    def arc_job_config(self, config, job_num, branches):
        config = law.ARCWorkflow.arc_job_config(self, config, job_num, branches)
        self._setup_render_variables(config, self.arc_workflow_requires())
        return config

    def arc_stageout_file(self):
        return law.util.rel_path(__file__, "files", "arc_stageout.sh")


class HTCondorWorkflow(law.htcondor.HTCondorWorkflow):
    """
    Batch systems are typically very heterogeneous by design, and so is HTCondor. Law does not aim
    to "magically" adapt to all possible HTCondor setups which would certainly end in a mess.
    Therefore we have to configure the base HTCondor workflow in law.contrib.htcondor to work with
    the VISPA environment. In most cases, like in this example, only a minimal amount of
    configuration is required.
    """
    htcondor_logs = luigi.BoolParameter()
    htcondor_gpus = luigi.IntParameter(default=2, significant=False, description="number "
        "of GPUs to request on the VISPA cluster")

    def htcondor_output_directory(self):
        # the directory where submission meta data should be stored
        return law.LocalDirectoryTarget(self.local_path())

    def htcondor_job_config(self, config, job_num, branches):
        # copy the entire environment
        config.custom_content.append(("getenv", "true"))
        # condor logs
        config.stdout = "out.txt"
        config.stderr = "err.txt"
        config.log = "log.txt"

        if self.htcondor_gpus > 0:
            config.custom_content.append(("request_gpus", self.htcondor_gpus))

        config.custom_content.append(("RequestMemory", "16000"))
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
        for subsystem in ["JetTaggingSF"]:
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


class AnalysisSandboxTask(law.SandboxTask):

    allow_empty_sandbox = True

    def sandbox_setup_cmds(self):
        cmds = super(AnalysisSandboxTask, self).sandbox_setup_cmds()

        cmds.append('export JTSF_CMSSW_SETUP="{}"'.format(os.environ["JTSF_CMSSW_SETUP"]))
        cmds.append("source {}".format(os.path.join(os.environ["JTSF_BASE"], "setup.sh")))
        cmds.append("source {}".format(os.path.join(
            os.environ["JTSF_BASE"], "singularity", "setup_$JTSF_DIST_VERSION.sh"))
        )

        return cmds

    def __init__(self, *args, **kwargs):
        super(AnalysisSandboxTask, self).__init__(*args, **kwargs)
        if self.sandbox_inst is not None:
            self.sandbox_inst.forward_env = False
            self.sandbox_inst.allow_binds = False


class UploadCMSSW(AnalysisTask, law.tasks.TransferLocalFile, AnalysisSandboxTask,
        law.cms.BundleCMSSW):

    force_upload = luigi.BoolParameter(default=False, description="force uploading")

    # settings for TransferLocalFile
    source_path = None

    version = None
    task_namespace = None

    def __init__(self, *args, **kwargs):
        super(UploadCMSSW, self).__init__(*args, **kwargs)

        self.has_run = False

    def store_parts(self):
        sl_dist_version = self.env["JTSF_DIST_VERSION"]
        return super(UploadCMSSW, self).store_parts() + (sl_dist_version,)

    def get_cmssw_path(self):
        return self.env["CMSSW_BASE"]

    def complete(self):
        if self.force_upload and not self.has_run:
            return False
        else:
            return super(UploadCMSSW, self).complete()

    def single_output(self):
        path = "{}.tgz".format(os.path.basename(self.get_cmssw_path()))
        return self.wlcg_target(path, fs="wlcg_fs_software")

    def output(self):
        return law.tasks.TransferLocalFile.output(self)

    def run(self):
        bundle = law.LocalFileTarget(is_tmp="tgz")
        self.bundle(bundle)
        self.transfer(bundle)

        self.has_run = True


class UploadSoftware(AnalysisTask, law.tasks.TransferLocalFile, AnalysisSandboxTask):

    version = None

    def __init__(self, *args, **kwargs):
        super(UploadSoftware, self).__init__(*args, **kwargs)
        self.source_path = self.env["JTSF_SOFTWARE"] + ".tgz"

    def store_parts(self):
        sl_dist_version = self.env["JTSF_DIST_VERSION"]
        return super(UploadSoftware, self).store_parts() + (sl_dist_version,)

    def single_output(self):
        return self.wlcg_target("software.tgz", fs="wlcg_fs_software")

    def run(self):
        # create the local bundle
        bundle = law.LocalFileTarget(self.source_path, is_tmp=True)
        def _filter(tarinfo):
            return None if re.search("(\.pyc|\/\.git|\.tgz)$", tarinfo.name) else tarinfo
        bundle.dump(os.path.splitext(self.source_path)[0], filter=_filter)
        self.publish_message("bundled software archive")

        # super run will upload all files for us
        super(UploadSoftware, self).run()


class UploadRepo(AnalysisTask, law.git.BundleGitRepository, law.tasks.TransferLocalFile):

    # settings for BundleGitRepository
    repo_path = os.environ["JTSF_BASE"]

    # settings for TransferLocalFile
    source_path = None

    version = None
    task_namespace = None

    def get_repo_path(self):
        return self.repo_path

    def single_output(self):
        path = "{}.{}.tgz".format(os.path.basename(self.get_repo_path()), self.checksum)
        return self.wlcg_target(path, fs="wlcg_fs_software")

    def output(self):
        return law.tasks.TransferLocalFile.output(self)

    def run(self):
        bundle = law.LocalFileTarget(is_tmp="tgz")
        self.bundle(bundle)
        self.transfer(bundle)
