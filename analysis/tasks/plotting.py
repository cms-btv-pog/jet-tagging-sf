# -*- coding: utf-8 -*-

import os
import luigi
import tarfile

from collections import defaultdict, OrderedDict

import numpy as np
import array
import itertools

import law
from law.parameter import CSVParameter
from law.target.local import LocalDirectoryTarget

from analysis.config.jet_tagging_sf import jes_total_shifts
from analysis.root import ROOTPlot
from analysis.tasks.base import AnalysisTask, WrapperTask
from analysis.tasks.hists import MergeHistograms, MergeScaleFactorWeights
from analysis.tasks.measurement import MeasureScaleFactors, MeasureCScaleFactors, FitScaleFactors


class PlotTask(AnalysisTask):
    def rebin_hist(self, hist, region, binning_type=None, b_tagger=None, truncate=False):
        if b_tagger is None:
            b_tagger = self.b_tagger

        if binning_type is None:
            bin_edges = [hist.GetBinLowEdge(idx) for idx in range(1, hist.GetNbinsX() + 2)]
        else:
            binning = self.config_inst.get_aux("binning")
            btagger_cfg = self.config_inst.get_aux("btaggers")[b_tagger]

            bin_edges = binning[region][b_tagger][binning_type]

        if truncate:
            # truncate < 0 bin
            bin_edges[0] = -0.1

        bin_edges = array.array("d", bin_edges)
        n_bins = len(bin_edges) - 1
        hist_rebinned = hist.Rebin(n_bins, "rebinned_{}".format(hist.GetName()), bin_edges)

        if truncate:
            # because of the truncation, the first bin content is filled into the underflow bin, fix this
            hist_rebinned.SetBinContent(1, hist.GetBinContent(1))
            hist_rebinned.SetBinError(1, hist.GetBinError(1))
        return hist_rebinned

    def output(self):
        return self.local_target("plots.tgz")


class PlotVariable(PlotTask):
    b_tagger = MergeHistograms.b_tagger
    iteration = MergeHistograms.iteration

    final_it = MergeHistograms.final_it

    category_tag = luigi.Parameter(default="merged")
    variable = CSVParameter(default=["jet{i_probe_jet}_{b_tag_var}_{region}_nominal"],
        description="Variable to plot, or multiple variables that are filled into one histogram. "
        "{} accesses auxiliary category information.")
    mc_split = luigi.ChoiceParameter(choices=["process", "flavor"], default="process")
    normalize = luigi.BoolParameter(description="Normalize MC histogram to data histogram")
    truncate = luigi.BoolParameter(description="Truncate the bin below zero, to be used "
        "for b-tag variable plots.")
    rebin = luigi.BoolParameter(description="Rebin variable to 'measurement' binning, only "
        "for b-tag variable plots. Not usable with category-optimized binning.")

    compare_iteration = luigi.IntParameter(default=-1, description="Secondary iteration to compare to. Can"
        "not be the final iteration.")
    logarithmic = luigi.BoolParameter(description="Plot y axis with logarithmic scale.")
    draw_stacked = luigi.BoolParameter(description="Plot MC processes separated by *mc_split*, "
        "combined in a stack.")

    mc_key = "mc"
    data_key = "data"

    def requires(self):
        reqs = {"hists": OrderedDict()}

        if self.compare_iteration >= 0:
            if self.compare_iteration == self.iteration:
                raise ValueError("Trying to compare identical iterations {} and {}".format(
                    self.iteration, self.compare_iteration))
            reqs["hists"]["secondary"] = MergeHistograms.req(self, branch=0, iteration=self.compare_iteration, final_it=False,
                version=self.get_version(MergeHistograms), _prefer_cli=["version"])

        reqs["hists"]["primary"] = MergeHistograms.req(self, branch=0, version=self.get_version(MergeHistograms),
            _prefer_cli=["version"])

        if self.normalize:
            reqs["scale"] = MeasureScaleFactors.req(self, iteration=0,
                version=self.get_version(MeasureScaleFactors), _prefer_cli=["version"])

        return reqs

    def associate_process(self, process):
        # associate process either to data or monte carlo
        # returns *add_to_data*, *sign* (1. or -1.)
        if process.is_data:
            return True, 1.
        else:
            return False, 1.

    def run(self):
        def add_hist(hist, new_hist, sign=1.):
            if hist is None:
                hist = new_hist.Clone()
                hist.Scale(sign)
            else:
                hist.Add(new_hist, sign)
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
            if category.has_tag((self.category_tag, self.b_tagger), mode=all):
                categories.append(category)

        # create plot objects
        plot_dict = {}
        limits_x = np.linspace(0., 1., len(inp["hists"]) + 1)
        for category in categories:
            plot = ROOTPlot(category.name, category.name)
            plot.create_pads(n_pads_x=len(inp["hists"]), n_pads_y=2, limits_x=limits_x, limits_y=[0., 0.3, 1.0])
            plot_dict[category] = plot

        for target_idx, inp_target in enumerate(inp["hists"].values()):
            with inp_target.load("r") as input_file:
                for category in categories:
                    data_hist = None
                    mc_hists = defaultdict(lambda: None)

                    for leaf_cat, _, children in category.walk_categories():
                        # we are only interested in leaves
                        if children:
                            continue

                        flavor = leaf_cat.get_aux("flavor", None)
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
                            for variable in self.variable:
                                # create variable name from template
                                aux = leaf_cat.aux.copy()
                                aux["b_tag_var"] = self.config_inst.get_aux("btaggers")[self.b_tagger]["variable"]
                                variable = variable.format(**aux)

                                hist = process_dir.Get(variable)

                                binning_type = "measurement" if self.rebin else None
                                hist = self.rebin_hist(hist, region, binning_type=binning_type, truncate=self.truncate)

                                add_to_data, sign = self.associate_process(process)
                                if add_to_data:
                                    data_hist = add_hist(data_hist, hist, sign=sign)
                                else:
                                    if self.normalize:  # apply "trigger" sfs as part of the normalization
                                        hist.Scale(scales[channel.name][region])

                                    key = process.name if self.mc_split == "process" else flavor
                                    mc_hists[key] = add_hist(mc_hists[key], hist, sign=sign)

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
                    hist_maximum = max([mc_hist_sum.GetMaximum(), data_hist.GetMaximum()])
                    plot = plot_dict[category]
                    # data and mc histograms
                    plot.cd(target_idx, 1)
                    if self.draw_stacked:
                        plot.draw(mc_hists, stacked=True, stack_maximum=1.5*hist_maximum)
                    else:
                        plot.draw({self.mc_key: mc_hist_sum}, line_color=None)
                    plot.draw({self.data_key: data_hist})

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
                    plot.draw_as_graph(ratio_mcerr_hist, options="2")
                    plot.draw({"data/mc": ratio_hist})

        for category, plot in plot_dict.items():
            plot.save(os.path.join(local_tmp.path,
                "{}_{}.pdf".format(category.name, self.variable)),
                draw_legend=True, log_y=self.logarithmic,
                lumi=self.config_inst.get_aux("lumi").values()[0]/1000.)
            del plot

        with outp.localize("w") as tmp:
            with tarfile.open(tmp.path, "w:gz") as tar:
                for plot_file in os.listdir(local_tmp.path):
                    tar.add(os.path.join(local_tmp.path, plot_file), arcname=plot_file)


class PlotContaminationEstimation(PlotVariable):
    data_key = "Data - nonZJets"
    mc_key = "ZJets"

    def associate_process(self, process):
        if process.is_data:
            return True, 1.

        if process.has_parent_process("dy"):
            return False, 1.
        else: # subtract non-ZJets from data
            return True, -1.


class PlotScaleFactor(PlotTask):
    hist_name = "sf"

    shifts = CSVParameter(default=["ALL"], description="Systematic shifts to plot."
        "Specify shift bases, up/down variations are created from this.")
    fix_normalization = FitScaleFactors.fix_normalization
    norm_to_nominal = luigi.BoolParameter()
    is_c_flavour = luigi.BoolParameter()

    b_taggers = CSVParameter(default=["deepcsv"])
    iterations = CSVParameter(default=[0])
    versions = CSVParameter(default=[None], description="Scale factor versions to compare."
        "The same version is used for all required tasks.")

    def __init__(self, *args, **kwargs):
        super(PlotScaleFactor, self).__init__(*args, **kwargs)
        # identifiers used in file names
        self.shifts_identifier = "_".join(self.shifts)
        self.file_identifiers = list(self.b_taggers) + [self.shifts_identifier]
        if self.is_c_flavour:
            self.file_identifiers.append("c")
        if self.norm_to_nominal:
            self.file_identifiers.append("normed")

        # Check if multiple shifts are present and thus have to be combined (envelope)
        self.multiple_shifts = "ALL" in self.shifts or len(self.shifts) >= 2

        if "ALL" in self.shifts:
            if self.is_c_flavour:
                self.shifts = MeasureCScaleFactors.shifts
            else:
                # do not double count jes shift (total shift and shift by sources)
                self.shifts = [shift for shift in MeasureScaleFactors.shifts
                    if not shift in jes_total_shifts]
        else:
            self.shifts = ["{}_{}".format(shift, direction) for shift, direction
                in itertools.product(self.shifts, ["up", "down"])]

        if not self.is_c_flavour:
            # make sure the nominal histograms are processed first
            if "nominal" in self.shifts:
                self.shifts.remove("nominal")
            self.shifts.insert(0, "nominal")

        if len(self.shifts) != len(set(self.shifts)):
            raise Exception("Duplicate shift in {}".format(self.shifts))

    def requires(self):
        reqs = OrderedDict()
        measure_task = MeasureCScaleFactors if self.is_c_flavour else MeasureScaleFactors
        for config in itertools.product(self.b_taggers, self.iterations, self.versions):
            b_tagger, iteration, version = config

            reqs[config] = OrderedDict()
            for shift in self.shifts:
                reqs[config][shift] = {
                    "fit": FitScaleFactors.req(self, shift=shift, b_tagger=b_tagger, iteration=iteration,
                        version=version if version is not None else self.get_version(FitScaleFactors),
                        _prefer_cli=["version"]),
                    "hist": measure_task.req(self, shift=shift, b_tagger=b_tagger, iteration=iteration,
                        version=version if version is not None else self.get_version(measure_task),
                        _prefer_cli=["version"])
                    }
            if self.fix_normalization and not self.is_c_flavour:
                reqs[config]["norm"] = MergeScaleFactorWeights.req(self, normalize_cerrs=False,
                    b_tagger=b_tagger, iteration=iteration,
                    version=version if version is not None else self.get_version(MergeScaleFactorWeights),
                    _prefer_cli=["version"])
        return reqs

    def output(self):
        filename = "plots_{}.tgz".format("_".join(self.file_identifiers))
        return self.local_target(filename)

    def run(self):
        import ROOT

        ROOT.PyConfig.IgnoreCommandLineOptions = True
        ROOT.gROOT.SetBatch()

        inp = self.input()
        outp = self.output()

        local_tmp = LocalDirectoryTarget(is_tmp=True)
        local_tmp.touch()

        plots = {}

        if self.norm_to_nominal and self.shifts[0] != "nominal":
            raise KeyError("'norm_to_nominal' is set to true, but no nominal values found.")

        for config, config_input in inp.items():
            b_tagger, iteration, version = config

            config_ids = [b_tagger]
            if len(self.iterations) > 1:
                config_ids.append("iteration {}".format(iteration))
            if len(self.versions) > 1:
                config_ids.append("version {}".format(version))
            config_id = ", ".join(config_ids)

            # get scaling factors for normalization
            if self.fix_normalization:
                norm_factors = config_input.pop("norm").load()["nominal"]

            nominal_fit_hists = {}
            # combined errors for multiple shifts
            up_shifted_fit_hists = defaultdict(dict)
            down_shifted_fit_hists = defaultdict(dict)

            for shift_idx, (shift, inp_target) in enumerate(config_input.items()):
                with inp_target["fit"]["sf"].load("r") as fit_file, \
                    inp_target["hist"]["scale_factors"].load("r") as hist_file:
                    for category_key in fit_file.GetListOfKeys():
                        category_name = category_key.GetName()
                        if not self.config_inst.has_category(category_name):
                            raise KeyError("Unknown category {}".format(category_name))

                        category = self.config_inst.get_category(category_name)
                        pt_range = category.get_aux("pt")
                        eta_range = category.get_aux("eta")
                        region = category.get_aux("region")

                        # same category name for different b-taggers
                        if len(self.b_taggers) > 1:
                            plot_category = category.name.replace("__" + b_tagger, "")
                        else:
                            plot_category = category.name

                        fit_category_dir = fit_file.Get(category_name)
                        fit_hist = fit_category_dir.Get(self.hist_name)

                        if shift == "nominal":
                            hist_category_dir = hist_file.Get(category_name)
                            hist = hist_category_dir.Get(self.hist_name)
                            # truncate first bin
                            hist = self.rebin_hist(hist, region, b_tagger=b_tagger, truncate=True)

                            # normalize histogram if required
                            if self.fix_normalization:
                                hist.Scale(norm_factors[category_name])

                            # make sure histograms are not cleaned up when the file is closed
                            nominal_fit_hists[plot_category] = fit_hist.Clone()
                            nominal_fit_hists[plot_category].SetDirectory(0)

                        # for c-jets, there is no nominal histogram
                        # Instead, all nominal values are set to 1
                        if shift_idx == 0 and self.is_c_flavour:
                            nominal_fit_hist = fit_hist.Clone()
                            for bin_idx in range(1, nominal_fit_hist.GetNbinsX() + 1):
                                nominal_fit_hist.SetBinContent(bin_idx, 1.0)

                            nominal_fit_hist.SetDirectory(0)
                            nominal_fit_hists[plot_category] = nominal_fit_hist

                        if shift != "nominal" and self.multiple_shifts:
                            # collect all shifted fit histograms to build envelope later
                            sys, direction = shift.rsplit("_", 1)
                            if direction == "up":
                                up_shifted_fit_hists[plot_category][sys] = fit_hist.Clone()
                                up_shifted_fit_hists[plot_category][sys].SetDirectory(0)
                            elif direction == "down":
                                down_shifted_fit_hists[plot_category][sys] = fit_hist.Clone()
                                down_shifted_fit_hists[plot_category][sys].SetDirectory(0)
                            else:
                                raise ValueError("Unknown direction {}".format(direction))

                        if self.norm_to_nominal:
                            fit_hist.Divide(nominal_fit_hists[plot_category])

                        # get same category key for all b-taggers
                        if plot_category in plots:
                            plot = plots[plot_category]
                        else:
                            plot = ROOTPlot(category.name, category.name)
                            plot.create_pads()
                            plots[plot_category] = plot
                        plot.cd(0, 0)
                        fit_hist.GetXaxis().SetRangeUser(-.1, 1.0)
                        y_min = 0.6 if self.norm_to_nominal else 0.
                        y_max = 1.4 if self.norm_to_nominal else 2.
                        fit_hist.GetYaxis().SetRangeUser(y_min, y_max)

                        if len(self.b_taggers) == 1:
                            title = self.config_inst.get_aux("btaggers")[b_tagger]["label"]
                        else:
                            title = "B-Tag Discriminant"
                        fit_hist.GetXaxis().SetTitle(title)

                        fit_hist.GetXaxis().SetTitleSize(.045)
                        fit_hist.GetYaxis().SetTitle("SF")
                        fit_hist.GetYaxis().SetTitleSize(.045)

                        if shift_idx == 0:
                            if not self.multiple_shifts or shift == "nominal":
                                # only draw this fit histogram if it is not part of a shifted envelope
                                plot.draw({"sf": fit_hist}, line_color=1, add_to_legend=False)
                            line = ROOT.TLine(0., 0., 0., 2.)
                            line.SetLineStyle(9)
                            plot.draw({"line": line}, add_same_option=False, line_color=1, add_to_legend=False)

                            # add category information to plot
                            if not np.isinf(pt_range[1]):
                                text = r"#splitline{%d < p_{T} < %d}{%.1f < |#eta| < %.1f}" % \
                                    (pt_range[0], pt_range[1], eta_range[0], eta_range[1])
                            else:
                                text = r"#splitline{p_{T} > %d}{%.1f < |#eta| < %.1f}" % \
                                    (pt_range[0], eta_range[0], eta_range[1])
                            plot.draw_text(text)
                        elif not self.multiple_shifts:
                            plot.draw({shift: fit_hist}, line_color=None)

                        if shift == "nominal" and not self.norm_to_nominal:
                            plot.draw({config_id + ", nominal": hist}, line_color=1)

            if self.multiple_shifts:
                for plot_category in plots:
                    plot = plots[plot_category]
                    plot.cd(0, 0)
                    errors_up = []
                    errors_down = []
                    for shift_idx, shift in enumerate(up_shifted_fit_hists[plot_category]):
                        up_shifted_hist = up_shifted_fit_hists[plot_category][shift]
                        down_shifted_hist = down_shifted_fit_hists[plot_category][shift]

                        for bin_idx in range(1, up_shifted_hist.GetNbinsX() + 1):
                            nominal_value = nominal_fit_hists[plot_category].GetBinContent(bin_idx)

                            # combine all shifts that have an effect in the same direction
                            # effect from <shift>_up/done systematics
                            diff_up = up_shifted_hist.GetBinContent(bin_idx) - nominal_value
                            diff_down = down_shifted_hist.GetBinContent(bin_idx) - nominal_value

                            # shift with effect in up/down direction
                            error_up = max([diff_up, diff_down, 0])
                            error_down = min([diff_up, diff_down, 0])

                            # detect systematics where up/down shift direction is the same
                            #if diff_up * diff_down > 0:
                            #    print "One sided shift: {}, {}".format(shift, plot_category)

                            # add in quadrature
                            if shift_idx == 0:
                                errors_up.append(error_up**2)
                                errors_down.append(error_down**2)
                            else:
                                errors_up[bin_idx - 1] += error_up**2
                                errors_down[bin_idx - 1] += error_down**2
                    errors_up = np.sqrt(errors_up)
                    errors_down = np.sqrt(errors_down)

                    # build shifted histograms
                    fit_hist_up = nominal_fit_hists[plot_category].Clone()
                    fit_hist_down = nominal_fit_hists[plot_category].Clone()

                    for bin_idx in range(1, fit_hist_up.GetNbinsX() + 1):
                        fit_hist_up.SetBinContent(bin_idx, fit_hist_up.GetBinContent(bin_idx)
                            + errors_up[bin_idx - 1])
                        fit_hist_down.SetBinContent(bin_idx, fit_hist_down.GetBinContent(bin_idx)
                            - errors_down[bin_idx - 1])

                    if self.norm_to_nominal:
                        fit_hist_up.Divide(nominal_fit_hists[plot_category])
                        fit_hist_down.Divide(nominal_fit_hists[plot_category])

                    plot.draw({config_id + ", up": fit_hist_up}, line_color=None)
                    plot.draw({config_id + ", down": fit_hist_down}, line_color=None)

        # save plots
        for plot_category in plots:
            plot = plots[plot_category]
            plot.save(os.path.join(local_tmp.path, "{}_{}.pdf".format(plot_category,
                self.shifts_identifier)), draw_legend=True,
                lumi=self.config_inst.get_aux("lumi").values()[0]/1000.)
            del plot

        with outp.localize("w") as tmp:
            with tarfile.open(tmp.path, "w:gz") as tar:
                for plot_file in os.listdir(local_tmp.path):
                    tar.add(os.path.join(local_tmp.path, plot_file), arcname=plot_file)


class PlotShiftedScaleFactorWrapper(AnalysisTask, law.WrapperTask):
    shifts = law.CSVParameter(default=[], description="shifts to require")
    skip_shifts = law.CSVParameter(default=[], description="shifts to skip, supports patterns")

    wrapped_task = PlotScaleFactor

    def __init__(self, *args, **kwargs):
        super(PlotShiftedScaleFactorWrapper, self).__init__(*args, **kwargs)

        if not self.shifts:
            # strip up/down direction from shifts, remove nominal variation
            shifts = MeasureScaleFactors.shifts
            self.shifts = list(set([shift.rsplit("_", 1)[0] for shift in shifts if not shift == "nominal"]))
        if self.skip_shifts:
            filter_fn = lambda d: not law.util.multi_match(d, self.skip_shifts)
            self.shifts = filter(filter_fn, self.shifts)

    def requires(self):
        def req(shift):
            return self.wrapped_task.req(self, shifts=[shift])

        return OrderedDict([(shift, req(shift)) for shift in self.shifts])
