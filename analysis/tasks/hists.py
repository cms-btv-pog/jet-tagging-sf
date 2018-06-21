# -*- coding: utf-8 -*-


import os
import collections

import law
import luigi
import six
from order.util import join_root_selection

from analysis.tasks.base import AnalysisTask, DatasetTask, GridWorkflow
from analysis.tasks.trees import MergeTrees, MergeMetaData


class WriteHistograms(DatasetTask, GridWorkflow, law.LocalWorkflow):

    file_merging = "trees"

    def workflow_requires(self):
        reqs = super(WriteHistograms, self).workflow_requires()

        if not self.cancel_jobs and not self.cleanup_jobs:
            reqs["meta"] = MergeMetaData.req(self, version=self.get_version(MergeMetaData),
                _prefer_cli=["version"])
            if not self.pilot:
                reqs["tree"] = MergeTrees.req(self, cascade_tree=-1,
                    version=self.get_version(MergeTrees), _prefer_cli=["version"])

        return reqs

    def requires(self):
        return {
            "tree": MergeTrees.req(self, cascade_tree=self.branch, branch=0,
                version=self.get_version(MergeTrees), _prefer_cli=["version", "workflow"]),
            "meta": MergeMetaData.req(self, version=self.get_version(MergeMetaData),
                _prefer_cli=["version"]),
        }

    def output(self):
        return self.wlcg_target("hists_{}.root".format(self.branch))

    def run(self):
        import ROOT

        inp = self.input()
        outp = self.output()
        outp.parent.touch(0o0770)

        # TODO:
        #  - weights, corrections, etc.
        #  - way to inject the tagging SF from previous iterations (requires some kind of parameter
        #    that is decremented in each iteration until e.g. 0)

        # get child categories
        categories = []
        channels = [self.config_inst.get_aux("dataset_channels")[self.dataset_inst]] \
            if self.dataset_inst.is_data else self.config_inst.channels.values()
        for channel in channels:
            for category, _, children in channel.walk_categories():
                if not children:
                    categories.append((channel, category))
        categories = list(set(categories))

        # get processes
        if len(self.dataset_inst.processes) != 1:
            raise NotImplementedError("cannot handle datasets with more than one process yet")
        processes = list(self.dataset_inst.processes.values())

        # build a progress callback
        progress = self.create_progress_callback(len(categories))

        # open the output file
        with outp.localize("w") as tmp:
            with tmp.dump("RECREATE") as output_file:
                with self.publish_step("creating root output file directories ..."):
                    process_dirs = {}
                    for _, category in categories:
                        output_file.cd()
                        category_dir = output_file.mkdir(category.name)
                        for process in processes:
                            category_dir.cd()
                            process_dir = category_dir.mkdir(process.name)
                            process_dir.Write()
                            process_dirs[(category.name, process.name)] = process_dir

                # open the input file and get the tree
                with inp["tree"].load("r") as input_file:
                    tree = input_file.Get("tree")
                    self.publish_message("{} events in tree".format(tree.GetEntries()))

                    # pt and eta aliases for jets and leptons
                    for obj in ["jet1", "jet2", "jet3", "jet4", "lep1", "lep2"]:
                        tree.SetAlias("{0}_pt".format(obj),
                            "({0}_px**2 + {0}_py**2)**0.5".format(obj))
                        tree.SetAlias("{0}_eta".format(obj),
                            "0.5 * log(({0}_E + {0}_pz) / ({0}_E - {0}_pz))".format(obj))

                    for i, (channel, category) in enumerate(categories):
                        self.publish_message("writing histograms in category {} ({}/{})".format(
                            category.name, i + 1, len(categories)))

                        for process in processes:
                            # weights
                            weights = []
                            if self.dataset_inst.is_mc:
                                # lumi weight
                                lumi = self.config_inst.get_aux("lumi")[channel]
                                x_sec = process.get_xsec(self.config_inst.campaign.ecm).nominal
                                sum_weights = inp["meta"].load()["event_weights"]["sum"]
                                lumi_weight = lumi * x_sec / sum_weights
                                weights.append(str(lumi_weight))

                            # totalWeight alias
                            while len(weights) < 2:
                                weights.insert(0, "1")
                            tree.SetAlias("totalWeight", join_root_selection(weights, op="*"))

                            # change into the correct directory
                            process_dirs[(category.name, process.name)].cd()

                            # actual projecting
                            for variable in self.config_inst.variables:
                                hist = ROOT.TH1F(variable.name, variable.full_title(root=True),
                                    *variable.binning)
                                hist.Sumw2()

                                # build the full selection string, including the total event weight
                                selection = category.selection
                                if variable.selection:
                                    selection = join_root_selection(selection, variable.selection)
                                selection = join_root_selection(selection, "{} > -10000".format(
                                    variable.expression))
                                selection = "({}) * totalWeight".format(selection)

                                # project and write the histogram
                                tree.Project(variable.name, variable.expression, selection)
                                hist.Write()

                        progress(i)


class MergeHistograms(AnalysisTask, law.CascadeMerge):

    merge_factor = 8  # TODO: optimize

    def create_branch_map(self):
        return law.CascadeMerge.create_branch_map(self)

    def cascade_workflow_requires(self):
        reqs = collections.OrderedDict()

        for dataset in self.config_inst.datasets:
            reqs[dataset.name] = WriteHistograms.req(self, dataset=dataset.name,
                version=self.get_version(WriteHistograms), _prefer_cli=["version", "workflow"])

        return reqs

    def trace_cascade_workflow_inputs(self, inputs):
        targets = []
        for d in inputs.values():
            targets.extend(list(d["collection"].targets.values()))
        return targets

    def cascade_requires(self, start_leaf, end_leaf):
        slices = []
        for dataset in self.config_inst.datasets:
            file_merging = WriteHistograms.file_merging
            n_files = self.config_inst.get_aux("get_file_merging")(file_merging, dataset)
            slices.extend([(dataset, i) for i in range(n_files)])

        target_slices = slices[start_leaf:end_leaf]
        return [
            WriteHistograms.req(self, dataset=d.name, branch=l,
                version=self.get_version(WriteHistograms), _prefer_cli=["version"])
            for d, l in target_slices
        ]

    def cascade_output(self):
        return self.wlcg_target("hists.root")

    def merge(self, inputs, output):
        tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
        tmp_dir.touch()

        with output.localize("w") as tmp_out:
            # fetch inputs
            with self.publish_step("fetching inputs ..."):
                def fetch(inp):
                    inp.copy_to_local(tmp_dir, cache=False)
                    return inp.basename

                def callback(i):
                    self.publish_message("fetch file {} / {}".format(i + 1, len(inputs)))

                bases = law.util.map_verbose(fetch, inputs, every=5, callback=callback)

            with self.publish_step("merging ..."):
                if len(bases) == 1:
                    self.publish_message("only 1 file to merge")
                    tmp_out.path = tmp_dir.child(bases[0]).path
                else:
                    # merge using hadd
                    bases = " ".join(bases)
                    cmd = "hadd -n 0 -d {} {} {}".format(tmp_dir.path, tmp_out.path, bases)
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
