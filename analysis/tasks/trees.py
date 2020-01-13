# -*- coding: utf-8 -*-


import os
import re
import sys
import json
import math
import shutil
import subprocess
import collections
import random

import law
import luigi
import six

from analysis.tasks.base import AnalysisTask, DatasetTask, WrapperTask, GridWorkflow, HTCondorWorkflow
from analysis.tasks.external import GetDatasetLFNs, DownloadSetupFiles
from analysis.util import wget, determine_xrd_redirector
from analysis.config.jet_tagging_sf import xrd_redirectors

class WriteTrees(DatasetTask, GridWorkflow, law.LocalWorkflow, HTCondorWorkflow):

    max_events = luigi.IntParameter(default=law.NO_INT)

    workflow_run_decorators = [law.decorator.notify]

    stream_input_file = False

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

        # create the temporary dir to run in
        tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
        tmp_dir.touch()

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

        # determine the xrd redirector and download the file
        if os.environ["JTSF_DIST_VERSION"] == "slc6":
            redirector = determine_xrd_redirector(lfn)
        else:
            redirector = random.choice(xrd_redirectors)
        xrd_url = "root://{}/{}".format(redirector, lfn)
        if self.stream_input_file:
            input_file = xrd_url
        else:
            dirname, basename = os.path.split(xrd_url)
            input_target = law.wlcg.WLCGFileTarget(basename,
                fs=law.wlcg.WLCGFileSystem(base=dirname, cache_config=None))
            local_input_target = tmp_dir.child("input_file.root", type="f")
            with self.publish_step("download input file ...", runtime=True):
                input_file = input_target.copy_to_local(local_input_target)
                input_file = law.target.file.add_scheme(input_file, "file")
                self.publish_message("size is {:.2f} {}".format(
                    *law.util.human_bytes(local_input_target.stat.st_size)))

        # cmsRun argument helper
        def cmsRunArg(key, value):
            return " ".join("{}={}".format(key, v) for v in law.util.make_list(value))

        output = self.output()
        # get global tag from dataset if defined, otherwise take default from config
        global_tag = self.dataset_inst.get_aux("global_tag", self.config_inst.get_aux("global_tag")[data_src])
        with output["tree"].localize("w") as tmp_tree, output["meta"].localize("w") as tmp_meta:
            args = [
                ("inputFiles", input_file),
                ("outputFile", tmp_tree.path),
                ("campaign", self.config_inst.campaign.name),
                ("metaDataFile", tmp_meta.path),
                ("isData", self.dataset_inst.is_data),
                ("globalTag", global_tag),
                ("lumiFile", setup_files["lumi_file"]),
                ("metFilters", self.config_inst.get_aux("metFilters")[data_src]),
                ("jesFiles", jes_files),
                ("jesRanges", jes_ranges),
                ("jesUncFiles", jes_unc_files),
                ("jesUncSrcFile", jes_unc_src_file),
                ("jesUncSources", self.config_inst.get_aux("jes_sources_{}".format(
                    self.config_inst.get_aux("jes_scheme"))
                )),
                ("jerPtResolutionFile", setup_files["jer_files"]["PtResolution"]),
                ("jerScaleFactorFile", setup_files["jer_files"]["SF"]),
                ("deepCSVWP", self.config_inst.get_aux("working_points")["deepcsv"]["medium"]),
                ("deepJetWP", self.config_inst.get_aux("working_points")["deepjet"]["medium"]),
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
            cfg_file = "treeMaker_cfg.py"
            cmd = "cmsRun " + law.util.rel_path(__file__, "files", cfg_file)
            cmd += " " + " ".join(cmsRunArg(*tpl) for tpl in args)

            # create environment
            env = os.environ.copy()
            env["CMSSW_SEARCH_PATH"] += ":" + setup_files_dir.path
            print("running command: {}".format(cmd))
            for obj in law.util.readable_popen(cmd, shell=True, executable="/bin/bash",
                    cwd=tmp_dir.path, env=env):
                if isinstance(obj, six.string_types):
                    print(obj)
                    if obj.startswith("Begin processing the"):
                        self._publish_message(obj)
                else:
                    if obj.returncode != 0:
                        raise Exception("cmsRun failed")

            if not tmp_tree.exists():
                raise Exception("output file not exising after cmsRun")


class WriteTreesWrapper(WrapperTask):

    wrapped_task = WriteTrees


class MergeTrees(DatasetTask, law.tasks.CascadeMerge, GridWorkflow, HTCondorWorkflow):

    merge_factor = 25

    workflow_run_decorators = [law.decorator.notify]

    def create_branch_map(self):
        return law.tasks.CascadeMerge.create_branch_map(self)

    def cascade_workflow_requires(self, **kwargs):
        return WriteTrees.req(self, version=self.get_version(WriteTrees), _prefer_cli=["version"],
            _exclude=("start_branch", "end_branch"), **kwargs)

    def cascade_requires(self, start_leaf, end_leaf):
        return [self.cascade_workflow_requires(branch=l) for l in range(start_leaf, end_leaf)]

    def trace_cascade_inputs(self, inputs):
        return [inp["tree"] for inp in inputs]

    def cascade_output(self):
        n_trees = self.config_inst.get_aux("get_file_merging")(self.config_inst,
            "trees", self.dataset)
        return law.SiblingFileCollection([
            self.wlcg_target("tree_{}.root".format(i)) for i in range(n_trees)
        ])

    def merge(self, inputs, output):
        tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
        tmp_dir.touch()

        # fetch inputs
        with self.publish_step("fetching inputs ..."):
            def fetch(inp):
                inp.copy_to_local(tmp_dir.child(inp.unique_basename, type="f"), cache=False)
                return inp.unique_basename

            def callback(i):
                self.publish_message("fetch file {} / {}".format(i + 1, len(inputs)))

            bases = law.util.map_verbose(fetch, inputs, every=5, callback=callback)

        # merge using hadd
        with self.publish_step("merging ..."):
            with output.localize("w", cache=False) as tmp_out:
                cmd = "hadd -O -n 0 -d {} {} {}".format(tmp_dir.path, tmp_out.path, " ".join(bases))
                code = law.util.interruptable_popen(cmd, shell="True", executable="/bin/bash",
                    cwd=tmp_dir.path)[0]
                if code != 0:
                    raise Exception("hadd failed")

                self.publish_message("merged file size: {:.2f} {}".format(
                    *law.util.human_bytes(os.stat(tmp_out.path).st_size)))

    def glite_output_postfix(self):
        return "_{}_{}".format(self.cascade_tree, self.cascade_depth)

    def arc_output_postfix(self):
        return self.glite_output_postfix()


class MergeTreesWrapper(WrapperTask):

    wrapped_task = MergeTrees

    cascade_tree = luigi.IntParameter(default=-1)


class MergeMetaData(DatasetTask):

    def requires(self):
        return WriteTrees.req(self, version=self.get_version(WriteTrees), _prefer_cli=["version"])

    def output(self):
        return self.wlcg_target("stats.json")

    @law.decorator.notify
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


class MergeMetaDataWrapper(WrapperTask):

    wrapped_task = MergeMetaData


class MeasureTreeSizes(AnalysisTask):

    merged_size = luigi.FloatParameter(default=1.0, description="target size of merged tree files "
        "in GB")

    def __init__(self, *args, **kwargs):
        super(MeasureTreeSizes, self).__init__(*args, **kwargs)

        self.has_run = False

    def complete(self):
        return self.has_run

    @law.decorator.notify
    def run(self):
        merged_files = collections.OrderedDict()

        total_size = 0.
        total_files = 0
        total_merged_files = 0

        cfg = law.config.Config.instance()

        for dataset in self.config_inst.datasets:
            print(" dataset {} ".format(dataset.name).center(80, "-"))

            # determine the full url to the remote directory, split it into uberftp door and path
            task = WriteTrees.req(self, dataset=dataset.name)
            url_path = task.output()["collection"].dir.uri()
            url_door = cfg.get("wlcg_fs", "base_listdir")
            m_path = re.match(r"^.+//(.*)\:\d+/.+(/pnfs/.+)$", url_path)
            m_door = re.match(r"^.+//(.*).+(/pnfs/.+)$", url_door)
            if (not m_path or not m_door):
                print("cannot parse url for dataset {}: {}, {}".format(
                    dataset.name, url_path, url_door)
                )
                continue
            _, path = m_path.groups()
            door, _ = m_door.groups()

            # run uberftp to fetch the sizes per file in bytes
            cmd = "uberftp -glob on {} 'ls {}/tree_*.root'".format(door, path)
            cmd += " | grep root | awk -F ' ' '{print $4}'"
            code, out, _ = law.util.interruptable_popen(cmd, shell=True, executable="/bin/bash",
                stdout=subprocess.PIPE)
            if code != 0:
                print("uberftp command failed")
                continue
            sizes = [int(s) for s in out.strip().split()]

            # calculate the number of files after merging
            n = len(sizes)
            sum_sizes = sum(sizes)
            if sum_sizes == 0:
                mean_size = -1
                merge_factor = -1
                merged_files[dataset.name] = -1
            else:
                mean_size = sum_sizes / float(n)
                target_size = self.merged_size * 1024.**3
                merge_factor = n if mean_size == 0 else min(n, int(round(target_size / mean_size)))
                merged_files[dataset.name] = int(math.ceil(n / float(merge_factor)))

            total_files += n
            total_size += sum_sizes
            total_merged_files += merged_files[dataset.name]

            print("files       : {} / {}".format(n, dataset.n_files))
            print("sum size    : {:.2f} {}".format(*law.util.human_bytes(sum_sizes)))
            print("mean size   : {:.2f} {}".format(*law.util.human_bytes(mean_size)))
            print("merge factor: {}".format(merge_factor))
            print("merged files: {}".format(merged_files[dataset.name]))

        # some output
        print(" summary ".center(80, "="))
        self.publish_message("total size   : {:.2f} {}".format(*law.util.human_bytes(total_size)))
        self.publish_message("total files  : {}".format(total_files))
        self.publish_message("after merging: {}".format(total_merged_files))
        self.publish_message("\nmerged files per dataset:")
        self.publish_message(json.dumps(merged_files, indent=4, separators=(",", ": ")))

        self.has_run = True
