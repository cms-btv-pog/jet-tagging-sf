# -*- coding: utf-8 -*-

import os
import subprocess
import collections
import shutil
import random

import law
import luigi
import six

from analysis.base import AnalysisTask, DatasetTask, GridWorkflow
from analysis.util import wget


class GetDatasetLFNs(DatasetTask, law.TransferLocalFile):

    source_path = None
    version = None

    def single_output(self):
        h = law.util.create_hash(sorted(self.dataset_inst.keys))
        return self.wlcg_target("lfns_{}.json".format(h))

    def run(self):
        lfns = []
        for key in self.dataset_inst.keys:
            print("get lfns for key {}".format(key))
            cmd = "dasgoclient -query='file dataset={}' -limit=0".format(key)
            code, out, _ = law.util.interruptable_popen(cmd, shell=True, stdout=subprocess.PIPE,
                executable="/bin/bash")
            if code != 0:
                raise Exception("dasgoclient query failed")
            lfns.extend(out.strip().split("\n"))

        tmp = law.LocalFileTarget(is_tmp="json")
        tmp.dump(lfns)
        self.transfer(tmp)


class DownloadSetupFiles(AnalysisTask, law.TransferLocalFile):

    source_path = None
    version = None

    def __init__(self, *args, **kwargs):
        super(DownloadSetupFiles, self).__init__(*args, **kwargs)

        self.source_files = self.get_source_files()

    def single_output(self):
        return self.wlcg_target("{}.tgz".format(law.util.create_hash(self.source_files)))

    @classmethod
    def create_path_hash(cls, path):
        if not isinstance(path, six.string_types):
            return None
        elif not path.startswith("http") and not path.startswith("/"):
            return None
        else:
            return "{}_{}".format(law.util.create_hash(path), os.path.basename(path))

    def get_source_files(self):
        files = {
            "lumi_file": self.config_inst.get_aux("lumi_file"),
            "normtag_file": self.config_inst.get_aux("normtag_file"),
        }
        return files

    def run(self):
        # create a tmp dir
        tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
        tmp_dir.touch()

        # download all setup files
        def download(src):
            h = self.create_path_hash(src)
            if h is None:
                return
            self.publish_message("download {}".format(src))
            dst = os.path.join(tmp_dir.path, h)
            if src.startswith("http"):
                wget(src, dst)
            else:
                shutil.copy2(src, dst)

        law.util.map_struct(download, self.source_files)

        # create a tmp archive
        tmp_arc = law.LocalFileTarget(is_tmp="tgz")
        tmp_arc.dump(tmp_dir)

        # transfer
        self.transfer(tmp_arc)

    def localize(self, **kwargs):
        # load the archive and unpack it into a temporary directory
        tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
        random.choice(self.output().targets).load(tmp_dir, **kwargs)

        def abspath(path):
            h = self.create_path_hash(path)
            return h and os.path.join(tmp_dir.path, h)

        return tmp_dir, law.util.map_struct(abspath, self.source_files)


class CreateTuples(DatasetTask, GridWorkflow, law.LocalWorkflow):

    def workflow_requires(self):
        return self.requires_from_branch()

    def requires(self):
        return {
            "lfns": GetDatasetLFNs.req(self, replicas=10),
            "files": DownloadSetupFiles.req(self, replicas=10),
        }

    def output(self):
        return self.wlcg_target("tuple_{}.root".format(self.branch))

    def run(self):
        lfn = random.choice(self.input()["lfns"].targets).load()[self.branch_data]
        setup_files_dir, setup_files = self.requires()["files"].localize()

        # gather information
        global_tag = self.config_inst.get_aux("global_tag")[self.dataset_inst.data_source]
        lumi_file = setup_files["lumi_file"]

        with self.output().localize("w") as tmp_output:
            # build the cmsRun command
            cmd = "cmsRun " + law.util.rel_path(__file__, "csvTreeMaker_cfg.py")
            cmd += " inputFiles=root://xrootd-cms.infn.it/{}".format(lfn)
            cmd += " outputFile={}".format(tmp_output.path)
            cmd += " isData={}".format(self.dataset_inst.is_data)
            cmd += " globalTag={}".format(global_tag)
            cmd += " lumiFile={}".format(lumi_file)
            cmd += " maxEvents=100"

            tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
            tmp_dir.touch()

            print("running command: {}".format(cmd))
            code = law.util.interruptable_popen(cmd, shell=True, executable="/bin/bash",
                cwd=tmp_dir.path)[0]
            if code != 0:
                raise Exception("cmsRun failed")

            if not tmp_output.exists():
                raise Exception("output file not exising after cmsRun")
