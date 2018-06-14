# -*- coding: utf-8 -*-

import law
import luigi
import six
from order.util import join_root_selection

from analysis.tasks.base import AnalysisTask, DatasetTask, GridWorkflow
from analysis.tasks.trees import MergeTrees, MergeMetaData


class WriteHistograms(DatasetTask, GridWorkflow, law.LocalWorkflow):

    def workflow_requires(self):
        if self.cancel_jobs or self.cleanup_jobs:
            return {}

        return self.requires_from_branch()

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

        # get processes
        if len(self.dataset_inst.processes) != 1:
            raise NotImplementedError("cannot handle datasets with more than one process yet")
        processes = list(self.dataset_inst.processes.values())

        # build a progress callback
        progress = self.create_progress_callback(len(categories))

        # open the output file
        with outp.dump("RECREATE") as output_file:
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
                        totalWeight = "1 * 1" if not weights else " * ".join(weights)
                        tree.SetAlias("totalWeight", totalWeight)

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
                                selection = join_root_selection(selection, category.selection)
                            selection = join_root_selection(selection, "{} > -10000".format(
                                variable.expression))
                            selection = "({}) * totalWeight".format(selection)

                            # project and write the histogram
                            tree.Project(variable.name, variable.expression, selection)
                            hist.Write()

                    progress(i)
