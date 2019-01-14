# -*- coding: utf-8 -*-

import os
import luigi
import tarfile

from collections import defaultdict, OrderedDict

import numpy as np
import array

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"figure.autolayout": True})
plt.style.use("ggplot")

from law.parameter import CSVParameter
from law.target.local import LocalDirectoryTarget

from analysis.root import ROOTPlot
from analysis.tasks.base import AnalysisTask
from analysis.tasks.hists import MergeHistograms
from analysis.tasks.measurement import MeasureScaleFactors, FitScaleFactors


class PlotTask(AnalysisTask):
    def rebin_hist(hist, region):
        # truncate < 0 bin
        binning = self.config_inst.get_aux("binning")
        btagger_cfg = self.config_inst.get_aux("btagger")

        bin_edges = array.array("d", binning[region][btagger_cfg["name"]]["measurement"])

        bin_edges[0] = -0.1
        n_bins = len(bin_edges) - 1
        hist_rebinned = hist.Rebin(n_bins, "rebinned_{}".format(hist.GetName()), bin_edges)
        # because of the truncation, the first bin content is filled into the underflow bin, fix this
        hist_rebinned.SetBinContent(1, hist.GetBinContent(1))
        hist_rebinned.SetBinError(1, hist.GetBinError(1))
        return hist_rebinned

    def output(self):
        return self.local_target("plots.tgz")


class PlotVariable(PlotTask):
    iteration = MergeHistograms.iteration
    final_it = MergeHistograms.final_it

    category_tag = luigi.Parameter(default="merged")
    variable = luigi.CSVParameter(default="jet{i_probe_jet}_deepcsv_bcomb_{region}_nominal",
        description="Variable to plot, or multiple variables that are filled into one histogram. "
        "{} accesses auxiliary category information.")
    mc_split = luigi.ChoiceParameter(choices=["process", "flavor"])
    normalize = luigi.BoolParameter()

    compare_it = luigi.IntParameter(default=-1, description="Secondary iteration to compare to. Can"
        "not be the final iteration.")

    def requires(self):
        reqs = {"hists": OrderedDict()}

        if self.compare_it >= 0:
            if self.compare_it == self.iteration:
                raise ValueError("Trying to compare identical iterations {} and {}".format(
                    self.iteration, self.compare_it))
            reqs["hists"]["secondary"] = MergeHistograms.req(self, branch=0, iteration=self.compare_it, final_it=False,
                version=self.get_version(MergeHistograms), _prefer_cli=["version"])

        reqs["hists"]["primary"] = MergeHistograms.req(self, branch=0, version=self.get_version(MergeHistograms),
            _prefer_cli=["version"])

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

        # create plot objects
        plot_dict = {}
        limits_x = np.linspace(0., 1., len(inp) + 1)
        for category in categories:
            plot = ROOTPlot(category.name, category.name)
            plot.create_pads(n_pads_x=len(inp), n_pads_y=2, limits_x=limits_x, limits_y=[0., 0.3, 1.0])
            plot_dict[category] = plot

        for target_idx, inp_target in enumerate(inp.values()):
            with inp_target.load("r") as input_file:
                for category in categories:
                    data_hist = None
                    mc_hists = defaultdict(lambda: None)

                    for leaf_cat, _, children in category.walk_categories():
                        # we are only interested in leaves
                        if children:
                            continue

                        flavor = leaf_cat.get_aux("flavor", None)
                        if self.normalize:
                            channel = leaf_cat.get_aux("channel")
                            region = leaf_cat.get_aux("region")

                        category_dir = input_file.GetDirectory(leaf_cat.name)
                        for process_key in category_dir.GetListOfKeys():
                            process = self.config_inst.get_process(process_key.GetName())
                            process_dir = category_dir.GetDirectory(process.name)

                            # avoid double counting of inclusive and flavor-dependent histograms
                            if flavor is not None:  # Not needed in case region isn't flavor specific
                                if process.is_data and flavor != "inclusive":
                                    continue
                                elif process.is_mc and flavor == "inclusive":
                                    continue
                            for variable in variables:
                                # create variable name from template
                                variable = self.variable.format(**leaf_cat.aux)

                                hist = process_dir.Get(variable)
                                if process.is_data:
                                    data_hist = add_hist(data_hist, hist)
                                else:
                                    if self.normalize:  # apply "trigger" sfs as aprt of the normalization
                                        hist.Scale(scales[channel.name][region])

                                    key = process.name if self.mc_split == "process" else flavor
                                    mc_hists[key] = add_hist(mc_hists[key], hist)

                    if self.normalize:  # normalize mc yield to data in this category
                        mc_yield = sum(hist.Integral() for hist in mc_hists.values())
                        data_yield = data_hist.Integral()
                        norm_factor = data_yield / mc_yield
                        for mc_hist in mc_hists.values():
                            mc_hist.Scale(norm_factor)

                    # get maximum value of hists/ stacks drawn to set axis ranges
                    mc_hist_sum = mc_hists.values()[0].Clone()
                    for mc_hist in mc_hists.values()[1:]:
                        mc_hist_sum.Add(mc_hist)
                    max_hist = mc_hist_sum.Clone() if \
                        (mc_hist_sum.GetMaximum() > data_hist.GetMaximum()) else data_hist.Clone()
                    max_hist.Scale(1.5)

                    plot = plot_dict[category]
                    # data and mc histograms
                    plot.cd(target_idx, 1)
                    plot.draw({"invis": max_hist}, invis=True)
                    plot.draw(mc_hists, stacked=True, options="SAME")
                    plot.draw({"data": data_hist}, options="SAME")

                    # ratio of data to mc below the main plot TODO: Error propagation
                    plot.cd(target_idx, 0)
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

        for category, plot in plot_dict.items():
            plot.save(os.path.join(local_tmp.path,
                "{}_{}.pdf".format(category.name, self.variable)))
            del plot

        with outp.localize("w") as tmp:
            with tarfile.open(tmp.path, "w:gz") as tar:
                for plot_file in os.listdir(local_tmp.path):
                    tar.add(os.path.join(local_tmp.path, plot_file), arcname=plot_file)


class PlotScaleFactor(AnalysisTask):
    hist_name = "sf"

    shifts = CSVParameter(default=["nominal"])
    iterations = CSVParameter(default=[0])
    fix_normalization = FitScaleFactors.fix_normalization

    def requires(self):
        reqs = OrderedDict()
        for iteration in self.iterations:
            reqs[iteration] = OrderedDict()
            for shift in self.shifts:
                reqs[iteration][shift] = {
                    "fit": FitScaleFactors.req(self, iteration=iteration, shift=shift,
                        version=self.get_version(FitScaleFactors), _prefer_cli=["version"]),
                    "hist": MeasureScaleFactors.req(self, iteration=iteration, shift=shift,
                        version=self.get_version(MeasureScaleFactors), _prefer_cli=["version"])
                    }
        return reqs

    def output(self):
        return self.local_target("plots.tgz")

    def run(self):
        import ROOT

        ROOT.PyConfig.IgnoreCommandLineOptions = True
        ROOT.gROOT.SetBatch()

        inp = self.input()
        outp = self.output()

        local_tmp = LocalDirectoryTarget(is_tmp=True)
        local_tmp.touch()

        plots = {}

        for iteration, inp_dict in inp.items():
            for shift, inp_target in inp_dict.items():
                with inp_target["fit"]["sf"].load("r") as fit_file, \
                    inp_target["hist"]["scale_factors"].load("r") as hist_file:
                    for category_key in fit_file.GetListOfKeys():
                        if not self.config_inst.has_category(category_key.GetName()):
                            continue

                        category = self.config_inst.get_category(category_key.GetName())
                        pt_range = category.get_aux("pt")
                        eta_range = category.get_aux("eta")
                        region = category.get_aux("region")
                        fit_category_dir = fit_file.Get(category_key.GetName())
                        hist_category_dir = hist_file.Get(category_key.GetName())

                        fit_hist = fit_category_dir.Get(self.hist_name)
                        hist = hist_category_dir.Get(self.hist_name)
                        hist = rebin_hist(hist, region)

                        if category in plots:
                            plot = plots[category]
                        else:
                            plot = ROOTPlot(category.name, category.name)
                            plot.create_pads(lumi=self.config_inst.get_aux("lumi").values()[0]/1000.)
                            plots[category] = plot
                        plot.cd(0, 0)
                        fit_hist.GetXaxis().SetRangeUser(-.1, 1.0)
                        fit_hist.GetYaxis().SetRangeUser(0., 2.0)
                        fit_hist.GetXaxis().SetTitle("DeepCSV")
                        fit_hist.GetXaxis().SetTitleSize(.045)
                        fit_hist.GetYaxis().SetTitle("SF")
                        fit_hist.GetYaxis().SetTitleSize(.045)

                        line = ROOT.TLine(0., 0., 0., 2.)
                        line.SetLineStyle(9)

                        plot.draw({"sf": fit_hist})
                        plot.draw({"hist": hist}, options="SAME")
                        plot.draw({"line": line})
                        # add category information to plot
                        if not np.isinf(pt_range[1]):
                            text = r"#splitline{%d < p_{T} < %d}{%.1f < |#eta| < %.1f}" % \
                                (pt_range[0], pt_range[1], eta_range[0], eta_range[1])
                        else:
                            text = r"#splitline{p_{T} > %d}{%.1f < |#eta| < %.1f}" % \
                                (pt_range[0], eta_range[0], eta_range[1])
                        plot.draw_text(text, .7, .8, size=0.04)
        # save plots
        for category in plots:
            plot = plots[category]
            plot.save(os.path.join(local_tmp.path, "{}.pdf".format(category.name)))
            del plot

        with outp.localize("w") as tmp:
            with tarfile.open(tmp.path, "w:gz") as tar:
                for plot_file in os.listdir(local_tmp.path):
                    tar.add(os.path.join(local_tmp.path, plot_file), arcname=plot_file)
