# -*- coding: utf-8 -*-

import os
import subprocess
import collections
import shutil
import random

import law
import luigi
import six

from analysis.tasks.base import AnalysisTask, DatasetTask, GridWorkflow
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
        h = law.util.create_hash(law.util.flatten(self.source_files))
        return self.wlcg_target("{}.tgz".format(h))

    @classmethod
    def create_path_hash(cls, path):
        if not isinstance(path, six.string_types):
            return None
        elif not path.startswith("http") and not path.startswith("/"):
            return None
        else:
            return "{}_{}".format(law.util.create_hash(path), os.path.basename(path))

    def get_source_files(self):
        # prepare JES files
        jes_url = lambda version, level: "https://raw.githubusercontent.com/cms-jet/JECDatabase" \
            "/master/textFiles/{0}/{0}_{1}_AK4PFchs.txt".format(version, level)
        jes_files = collections.defaultdict(lambda: collections.defaultdict(dict))
        for src in ("mc", "data"):
            for _, _, version in self.config_inst.get_aux("jes_version")[src]:
                for level in self.config_inst.get_aux("jes_levels")[src] + ["Uncertainty"]:
                    jes_files[src][version][level] = jes_url(version, level)
        jes_unc_src_file = jes_url(self.config_inst.get_aux("jes_version")["mc"][0][2],
            "UncertaintySources")

        # prepare JER files
        jer_url = lambda version, src: "https://raw.githubusercontent.com/cms-jet/JRDatabase" \
            "/master/textFiles/{0}/{0}_{1}_AK4PFchs.txt".format(version, src)
        jer_files = {}
        for src in ("SF", "PtResolution", "PhiResolution"):
            jer_files[src] = jer_url(self.config_inst.get_aux("jer_version") + "_MC", src)

        return {
            "lumi_file": self.config_inst.get_aux("lumi_file"),
            "normtag_file": self.config_inst.get_aux("normtag_file"),
            "jes_files": jes_files,
            "jes_unc_src_file": jes_unc_src_file,
            "jer_files": jer_files,
        }

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


class CreateTrees(DatasetTask, GridWorkflow, law.LocalWorkflow):

    def workflow_requires(self):
        return self.requires_from_branch()

    def requires(self):
        return {
            "lfns": GetDatasetLFNs.req(self, replicas=10),
            "files": DownloadSetupFiles.req(self, replicas=10),
        }

    def output(self):
        return self.wlcg_target("tree_{}.root".format(self.branch))

    def run(self):
        lfn = random.choice(self.input()["lfns"].targets).load()[self.branch_data]
        setup_files_dir, setup_files = self.requires()["files"].localize()

        # manage jes files
        data_src = self.dataset_inst.data_source
        jes_versions = self.config_inst.get_aux("jes_version")[data_src]
        jes_levels = self.config_inst.get_aux("jes_levels")[data_src]
        jes_ranges = law.util.flatten(tpl[:2] for tpl in jes_versions)
        jes_files_dict = setup_files["jes_files"][data_src]
        jes_files = law.util.flatten([
            [jes_files_dict[version][level] for level in jes_levels]
            for _, _, version in jes_versions
        ])
        jes_unc_files = [jes_files_dict[version]["Uncertainty"] for _, _, version in jes_versions]
        jes_unc_src_file = setup_files["jes_unc_src_file"] if self.dataset_inst.is_mc else ""

        # cmsRun argument helper
        def cmsRunArg(key, value):
            return " ".join("{}={}".format(key, v) for v in law.util.make_list(value))

        with self.output().localize("w") as tmp_output:
            args = [
                ("inputFiles", "root://xrootd-cms.infn.it/{}".format(lfn)),
                ("outputFile", tmp_output.path),
                ("isData", self.dataset_inst.is_data),
                ("globalTag", self.config_inst.get_aux("global_tag")[data_src]),
                ("lumiFile", setup_files["lumi_file"]),
                ("metFilters", self.config_inst.get_aux("metFilters")[data_src]),
                ("jesFiles", jes_files),
                ("jesRanges", jes_ranges),
                ("jesUncFiles", jes_unc_files),
                ("jesUncSrcFile", jes_unc_src_file),
                ("jesUncSources", self.config_inst.get_aux("jes_sources")),
                ("maxEvents", 20000)
            ]

            # triggers
            for channel_inst, triggers in self.config_inst.get_aux("triggers").items():
                # special rules may apply for real datasets as triggers can be run dependent
                if self.dataset_inst.is_data:
                    d_ch = self.config_inst.get_aux("dataset_channels")[self.dataset_inst]
                    if d_ch == channel_inst:
                        triggers = self.config_inst.get_aux("data_triggers").get(
                            self.dataset_inst, triggers)
                args.append((channel_inst.name + "Triggers", triggers))

            # lepton channel for data
            if self.dataset_inst.is_data:
                ch = self.config_inst.get_aux("dataset_channels")[self.dataset_inst].name
                args.append(("leptonChannel", ch))

            # build the cmsRun command
            cmd = "cmsRun " + law.util.rel_path(__file__, "../cmssw/treeMaker_ICHEP18_cfg.py")
            cmd += " " + " ".join(cmsRunArg(*tpl) for tpl in args)

            tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
            tmp_dir.touch()

            print("running command: {}".format(cmd))
            code = law.util.interruptable_popen(cmd, shell=True, executable="/bin/bash",
                cwd=tmp_dir.path)[0]
            if code != 0:
                raise Exception("cmsRun failed")

            if not tmp_output.exists():
                raise Exception("output file not exising after cmsRun")
