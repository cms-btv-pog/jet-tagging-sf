# -*- coding: utf-8 -*-

import os
import luigi
import tarfile

from collections import defaultdict

from law.parameter import CSVParameter
from law.target.local import LocalDirectoryTarget

from analysis.root import ROOTPlot
from analysis.tasks.base import AnalysisTask
from analysis.tasks.hists import MergeHistograms
from analysis.tasks.measurement import MeasureScaleFactors

class PlotTask(AnalysisTask):
    iteration = MergeHistograms.iteration

    def store_parts(self):
        return super(PlotTask, self).store_parts() + (self.iteration,)

    def output(self):
        return self.local_target("plots.tgz")


class PlotVariable(PlotTask):
    category_tag = luigi.Parameter(default="merged")
    variable = luigi.Parameter(default="jet{i_probe_jet}_deepcsv_bcomb")
    mc_split = luigi.ChoiceParameter(choices=["process", "flavor"])
    normalize = luigi.BoolParameter()

    def requires(self):
        reqs = {
            "hist": MergeHistograms.req(self, branch=0, version=self.get_version(MergeHistograms),
                _prefer_cli=["version"])
        }
        if self.normalize:
            reqs["scale"] = MeasureScaleFactors.req(self, iteration=0,
                version=self.get_version(MeasureScaleFactors), _prefer_cli=["version"])

        return reqs

    def run(self):
        def add_hist(hist, new_hist):
            if hist is None:
                hist = new_hist.Clone()
            else:
                hist.Add(new_hist)
            return hist

        import ROOT

        ROOT.PyConfig.IgnoreCommandLineOptions = True
        ROOT.gROOT.SetBatch()

        inp = self.input()
        outp = self.output()

        if self.normalize:
            scales = inp["scale"]["channel_scales"].load()

        local_tmp = LocalDirectoryTarget(is_tmp=True)
        local_tmp.touch()

        categories = []
        for category, _, _ in self.config_inst.walk_categories():
            if category.has_tag(self.category_tag):
                categories.append(category)

        with inp["hist"].load("r") as input_file:
            for category in categories:
                data_hist = None
                mc_hists = defaultdict(lambda: None)

                for leaf_cat, _, children in category.walk_categories():
                    # we are only interested in leaves
                    if children:
                        continue

                    # create variable name from template
                    variable = self.variable.format(**leaf_cat.aux)

                    flavor = leaf_cat.get_aux("flavor")
                    channel = leaf_cat.get_aux("channel")
                    region = leaf_cat.get_aux("region")

                    category_dir = input_file.GetDirectory(leaf_cat.name)
                    for process_key in category_dir.GetListOfKeys():
                        process = self.config_inst.get_process(process_key.GetName())
                        process_dir = category_dir.GetDirectory(process.name)

                        # avoid double counting of inclusive and flavor-dependent histograms
                        if process.is_data and flavor != "inclusive":
                            continue
                        elif process.is_mc and flavor == "inclusive":
                            continue

                        hist = process_dir.Get(variable)
                        if process.is_data:
                            data_hist = add_hist(data_hist, hist)
                        else:
                            if self.normalize:
                                hist.Scale(scales[channel.name][region])
                            key = process.name if self.mc_split == "process" else flavor
                            mc_hists[key] = add_hist(mc_hists[key], hist)

                # get maximum value of hists/ stacks drawn to set axis ranges
                mc_hist_sum = mc_hists.values()[0].Clone()
                for mc_hist in mc_hists.values()[1:]:
                    mc_hist_sum.Add(mc_hist)
                max_hist = mc_hist_sum.Clone() if \
                    (mc_hist_sum.GetMaximum() > data_hist.GetMaximum()) else data_hist.Clone()
                max_hist.Scale(1.5)

                # data and mc histograms
                plot = ROOTPlot(category.name, category.name)
                plot.create_pads(n_pads_y=2, limits_y=[0., 0.3, 1.0])
                plot.cd(0, 1)
                plot.draw({"invis": max_hist}, invis=True)
                plot.draw(mc_hists, stacked=True, options="SAME")
                plot.draw({"data": data_hist}, options="SAME")

                # ratio of data to mc below the main plot TODO: Error propagation
                plot.cd(0, 0)
                # mc error band
                ratio_mcerr_hist = mc_hist_sum.Clone()
                ratio_mcerr_hist.Divide(mc_hist_sum)
                # ratio
                ratio_hist = data_hist.Clone()
                ratio_hist.Divide(mc_hist_sum)
                ratio_hist.GetYaxis().SetRangeUser(0.5, 1.5)
                plot.draw({"invis": ratio_hist}, invis=True)
                plot.draw_as_graph(ratio_mcerr_hist, options=["SAME", "2"])
                plot.draw({"data/mc": ratio_hist}, options="SAME")

                plot.save(os.path.join(local_tmp.path,
                    "{}_{}.pdf".format(category.name, self.variable)))
                del plot

        with outp.localize("w") as tmp:
            with tarfile.open(tmp.path, "w:gz") as tar:
                for plot_file in os.listdir(local_tmp.path):
                    tar.add(os.path.join(local_tmp.path, plot_file), arcname=plot_file)


class PlotScaleFactor(PlotTask):

    hist_name = "sf"

    def requires(self):
        return MeasureScaleFactors.req(self, version=self.get_version(MeasureScaleFactors),
                _prefer_cli=["version"])

    def run(self):
        import ROOT

        ROOT.PyConfig.IgnoreCommandLineOptions = True
        ROOT.gROOT.SetBatch()

        inp = self.input()
        outp = self.output()

        local_tmp = LocalDirectoryTarget(is_tmp=True)
        local_tmp.touch()

        with inp["scale_factors"].load("r") as input_file:
            for category_key in input_file.GetListOfKeys():
                if not self.config_inst.has_category(category_key.GetName()):
                    continue

                category = self.config_inst.get_category(category_key.GetName())
                category_dir = input_file.Get(category_key.GetName())

                hist = category_dir.Get(self.hist_name)

                plot = ROOTPlot(category.name, category.name)
                plot.create_pads()
                plot.cd(0, 0)
                hist.GetYaxis().SetRangeUser(0., 2.0)
                plot.draw({"sf": hist})
                plot.save(os.path.join(local_tmp.path, "{}.pdf".format(category.name)))
                del plot

        with outp.localize("w") as tmp:
            with tarfile.open(tmp.path, "w:gz") as tar:
                for plot_file in os.listdir(local_tmp.path):
                    tar.add(os.path.join(local_tmp.path, plot_file), arcname=plot_file)
