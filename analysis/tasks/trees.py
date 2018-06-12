# -*- coding: utf-8 -*-

import os
import subprocess
import collections
import shutil

import law
import luigi
import six

from analysis.tasks.base import AnalysisTask, DatasetTask, GridWorkflow
from analysis.util import wget, determine_xrd_redirector


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
        self.output().random_target().load(tmp_dir, **kwargs)

        def abspath(path):
            h = self.create_path_hash(path)
            return h and os.path.join(tmp_dir.path, h)

        return tmp_dir, law.util.map_struct(abspath, self.source_files)


class WriteTrees(DatasetTask, GridWorkflow, law.LocalWorkflow):

    max_events = luigi.IntParameter(default=law.NO_INT)

    def workflow_requires(self):
        if self.cancel_jobs or self.cleanup_jobs:
            return {}

        return self.requires_from_branch()

    def requires(self):
        return {
            "lfns": GetDatasetLFNs.req(self, replicas=10),
            "files": DownloadSetupFiles.req(self, replicas=10),
        }

    def output(self):
        return {
            "tree": self.wlcg_target("tree_{}.root".format(self.branch)),
            "meta": self.wlcg_target("meta_{}.root".format(self.branch)),
        }

    def run(self):
        lfn = self.input()["lfns"].random_target().load()[self.branch_data]
        setup_files_dir, setup_files = self.requires()["files"].localize()

        # determine the xrd redirector
        redirector = determine_xrd_redirector(lfn)

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

        output = self.output()
        with output["tree"].localize("w") as tmp_tree, output["meta"].localize("w") as tmp_meta:
            args = [
                ("inputFiles", "root://{}/{}".format(redirector, lfn)),
                ("outputFile", tmp_tree.path),
                ("metaDataFile", tmp_meta.path),
                ("isData", self.dataset_inst.is_data),
                ("globalTag", self.config_inst.get_aux("global_tag")[data_src]),
                ("lumiFile", setup_files["lumi_file"]),
                ("metFilters", self.config_inst.get_aux("metFilters")[data_src]),
                ("jesFiles", jes_files),
                ("jesRanges", jes_ranges),
                ("jesUncFiles", jes_unc_files),
                ("jesUncSrcFile", jes_unc_src_file),
                ("jesUncSources", self.config_inst.get_aux("jes_sources")),
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

            # max events
            if not law.is_no_param(self.max_events):
                args.append(("maxEvents", self.max_events))

            # build the cmsRun command
            cfg_file = "treeMaker_{}_cfg.py".format(os.getenv("JTSF_CMSSW_SETUP"))
            cmd = "cmsRun " + law.util.rel_path(__file__, cfg_file)
            cmd += " " + " ".join(cmsRunArg(*tpl) for tpl in args)

            tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
            tmp_dir.touch()

            print("running command: {}".format(cmd))
            code = law.util.interruptable_popen(cmd, shell=True, executable="/bin/bash",
                cwd=tmp_dir.path)[0]
            if code != 0:
                raise Exception("cmsRun failed")

            if not tmp_tree.exists():
                raise Exception("output file not exising after cmsRun")


class MergeTrees(DatasetTask, law.CascadeMerge):

    merge_factor = 8

    def create_branch_map(self):
        return law.CascadeMerge.create_branch_map(self)

    def cascade_workflow_requires(self, **kwargs):
        return WriteTrees.req(self, version=self.get_version(WriteTrees), _prefer_cli=["version"],
            **kwargs)

    def cascade_requires(self, start_leaf, end_leaf):
        return [self.cascade_workflow_requires(branch=l) for l in range(start_leaf, end_leaf)]

    def trace_cascade_inputs(self, inputs):
        return [inp["tree"] for inp in inputs]

    def cascade_output(self):
        n_trees = self.config_inst.get_aux("file_merging")["trees"].get(self.dataset, 1)
        return law.SiblingFileCollection([
            self.wlcg_target("tree_{}.root".format(i)) for i in range(n_trees)
        ])

    def merge(self, inputs, output):
        tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
        tmp_dir.touch()

        # fetch inputs
        with self.publish_step("fetching inputs ..."):
            def fetch(inp):
                inp.copy_to_local(tmp_dir, cache=False)
                return inp.basename

            def callback(i):
                self.publish_message("fetch file {} / {}".format(i + 1, len(inputs)))

            bases = law.util.map_verbose(fetch, inputs, every=5, callback=callback)

        # merge using hadd
        with self.publish_step("merging ..."):
            with output.localize("w", cache=False) as tmp_out:
                cmd = "hadd -O -n 0 -d {} {} {}".format(tmp_dir.path, tmp_out.path, " ".join(bases))
                code, _, _ = law.util.interruptable_popen(cmd, shell="True", executable="/bin/bash",
                    cwd=tmp_dir.path)
                if code != 0:
                    raise Exception("hadd failed")

                self.publish_message("merged file size: {:.2f} {}".format(
                    *law.util.human_bytes(os.stat(tmp_out.path).st_size)))

    def glite_output_postfix(self):
        return "_{}_{}".format(self.cascade_tree, self.cascade_depth)

    def arc_output_postfix(self):
        return self.glite_output_postfix()


class MergeMetaData(DatasetTask):

    def requires(self):
        return WriteTrees.req(self, version=self.get_version(WriteTrees), _prefer_cli=["version"])

    def output(self):
        return self.wlcg_target("stats.json")

    def run(self):
        stats = {}

        def load(inp):
            with inp["meta"].load(formatter="root", cache=False) as tfile:
                for key in tfile.GetListOfKeys():
                    key = key.GetName()
                    obj = tfile.Get(key)
                    if not obj.__class__.__name__.startswith("TH1"):
                        continue

                    n_bins = obj.GetNbinsX()
                    values = stats.setdefault(key, n_bins * [0.])
                    for i in range(n_bins):
                        values[i] += obj.GetBinContent(i + 1)

        def callback(i):
            self.publish_message("loading meta data {} / {}".format(i + 1, len(coll)))
            self.publish_progress(100. * (i + 1) / len(coll))

        coll = self.input()["collection"]
        law.util.map_verbose(load, coll.targets.values(), every=25, callback=callback)

        # manual formatting of some entries
        for key in ["events", "event_weights", "selected_events", "selected_event_weights"]:
            stats[key] = {
                "sum": stats[key][0] + stats[key][1],
                "positive": stats[key][1],
                "negative": stats[key][0],
            }

        self.output().dump(stats, formatter="json", indent=4)
