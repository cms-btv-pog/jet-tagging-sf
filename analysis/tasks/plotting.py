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


class PlotVariable(AnalysisTask):
    category_tag = luigi.Parameter(default="merged")
    variables = CSVParameter(default=["jet1_pt"])

    def requires(self):
        return MergeHistograms.req(self, branch=0, version=self.get_version(MergeHistograms),
                _prefer_cli=["version"])

    def output(self):
        return self.local_target("plots.tgz")

    def run(self):
        def add_hist(hist, new_hist):
            if hist is None:
                hist = new_hist.Clone()
            else:
                hist.Add(new_hist)
            return hist

        import ROOT

        inp = self.input()
        outp = self.output()

        local_tmp = LocalDirectoryTarget(is_tmp=True)
        local_tmp.touch()

        categories = []
        for category, _, _ in self.config_inst.walk_categories():
            if category.has_tag(self.category_tag):
                categories.append(category)

        with inp.load("r") as input_file:
            for category in categories:

                for variable in self.variables:
                    data_hist = None
                    mc_hists = defaultdict(lambda: None)

                    for leaf_cat, _, children in category.walk_categories():
                        # we are only interested in leaves
                        if children:
                            continue
                        category_dir = input_file.GetDirectory(leaf_cat.name)

                        for process_key in category_dir.GetListOfKeys():
                            process = self.config_inst.get_process(process_key.GetName())
                            process_dir = category_dir.GetDirectory(process.name)

                            hist = process_dir.Get(variable)
                            if process.is_data:
                                data_hist = add_hist(data_hist, hist)
                            else:
                                mc_hists[process] = add_hist(mc_hists[process], hist)

                    plot = ROOTPlot(category.name, category.name)
                    plot.draw(mc_hists.values(), stacked=True)
                    plot.draw(data_hist, options="SAME")

                    plot.save(os.path.join(local_tmp.path,
                        "{}_{}.pdf".format(category.name, variable)))
                    del plot

        with outp.localize("w") as tmp:
            with tarfile.open(tmp.path, "w:gz") as tar:
                for plot_file in os.listdir(local_tmp.path):
                    tar.add(os.path.join(local_tmp.path, plot_file), arcname=plot_file)


class PlotScaleFactor(AnalysisTask):

    hist_name = "sf"

    def requires(self):
        return MeasureScaleFactors.req(self, version=self.get_version(MeasureScaleFactors),
                _prefer_cli=["version"])

    def output(self):
        return self.local_target("plots.tgz")

    def run(self):
        import ROOT

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

                plot = ROOTPlot()
                plot.draw(hist)
                plot.save(os.path.join(local_tmp.path, "{}.pdf".format(category.name)))

        with outp.localize("w") as tmp:
            with tarfile.open(tmp.path, "w:gz") as tar:
                for plot_file in os.listdir(local_tmp.path):
                    tar.add(os.path.join(local_tmp.path, plot_file), arcname=plot_file)
