# -*- coding: utf-8 -*-


import os
import array
import luigi
import itertools
import numpy as np

from collections import defaultdict

from analysis.config.jet_tagging_sf import jes_sources, jes_total_shifts
from analysis.tasks.base import AnalysisTask, ShiftTask, WrapperTask
from analysis.tasks.hists import MergeHistograms, GetScaleFactorWeights, MergeScaleFactorWeights
from analysis.tasks.util import OptimizeBinning
from analysis.util import format_shifts


class MeasureScaleFactors(ShiftTask):

    iteration = MergeHistograms.iteration
    b_tagger = MergeHistograms.b_tagger
    optimize_binning = MergeHistograms.optimize_binning
    category_tags = MergeHistograms.category_tags

    shifts = {"nominal"} | format_shifts(jes_sources, prefix="jes") | \
        format_shifts(["lf", "hf", "lf_stats1", "lf_stats2", "hf_stats1", "hf_stats2"])

    def requires(self):
        reqs = {
            "hist": MergeHistograms.req(self, branch=0, version=self.get_version(MergeHistograms),
                _prefer_cli=["version"])
        }
        if self.iteration > 0 or self.effective_shift != "nominal":
            reqs["scale"] = MeasureScaleFactors.req(self, iteration=0, shift="nominal",
                version=self.get_version(MeasureScaleFactors), _prefer_cli=["version"])
        if self.optimize_binning:
            reqs["binning"] = OptimizeBinning.req(self, version=self.get_version(OptimizeBinning),
                _prefer_cli=["version"])
        return reqs

    def store_parts(self):
        binning_part = "optimized" if self.optimize_binning else "default"
        return super(MeasureScaleFactors, self).store_parts() + (self.b_tagger,) + (self.iteration,) \
            + (binning_part,)

    def output(self):
        outputs = {"scale_factors": self.wlcg_target("scale_factors.root")}
        if self.iteration == 0 and self.shift == "nominal":
            outputs["channel_scales"] = self.wlcg_target("channel_scales.json")
        return outputs

    def get_flavor_component(self, flavor, region):
        """
        Depending on the region, c jets are counted for either the ``"heavy"`` or the ``"light"``
        flavor component.
        """
        if region == "hf":
            return "heavy" if flavor == "b" else "light"
        elif region == "lf":
            return "light" if flavor == "udsg" else "heavy"
        else:
            raise ValueError("unexpected region %s" % region)

    def run(self):
        def hist_integral(hist):
            return hist.Integral(0, hist.GetNbinsX() + 1)

        inp = self.input()
        outp = self.output()

        # get categories in which we measure the scale factors
        # these are stored in the config itself as we measure them inclusively over channels
        categories = []
        for category, _, _ in self.config_inst.walk_categories():
            if category.has_tag(("merged", self.b_tagger), mode=all) and category.get_aux("phase_space") == "measure":
                if len(self.category_tags) > 0 and not category.has_tag(self.category_tags, mode=any):
                    continue
                categories.append(category)

        # get categories from which to determine the rate scaling of MC to data
        # only needed for the first iteration, where the scaling is saved for further use
        if self.iteration == 0 and self.shift == "nominal":
            scale_categories = {}
            for channel in self.config_inst.channels:
                scale_categories[channel] = {}
                for category, _, children in channel.walk_categories():
                    if category.has_tag(("scales", self.b_tagger), mode=all) and category.get_aux("phase_space") == "closure":
                        if category.get_aux("config", None) != self.config_inst.name:
                            continue

                        region = category.get_aux("region")
                        scale_categories[channel][region] = category

        btagger_cfg = self.config_inst.get_aux("btaggers")[self.b_tagger]

        # get category-dependent binning if optimized binning is used
        if self.optimize_binning:
            category_binnings = inp["binning"].load()

        # category -> component (heavy/light) -> histogram
        hist_dict = {}
        # category -> histogram
        sf_dict = {}

        with inp["hist"].load("r") as input_file:
            # get scale factor to scale MC (withouts b-tag SFs) to data per channel
            if self.iteration == 0 and self.shift == "nominal":
                variable_name = "jet1_{}_bcomb".format(self.b_tagger)
                scales = defaultdict(dict)
                for channel, region_categories in scale_categories.items():
                    for region, category in region_categories.items():
                        data_yield, mc_yield = 0., 0.
                        category_dir = input_file.GetDirectory(category.name)

                        # sum over processes to get mc and data yields
                        for process_key in category_dir.GetListOfKeys():
                            process = self.config_inst.get_process(process_key.GetName())
                            process_dir = category_dir.GetDirectory(process.name)
                            hist = process_dir.Get("{}_{}_{}".format(variable_name, region, self.shift))
                            if process.is_data:
                                data_yield += hist_integral(hist)
                            else:
                                mc_yield += hist_integral(hist)
                        scale = data_yield / mc_yield
                        scales[channel.name][region] = scale
            else:
                scales = inp["scale"]["channel_scales"].load()
            for category in categories:
                region = category.get_aux("region")

                hist_dict[category] = {}
                for leaf_cat, _, children in category.walk_categories():
                    # we are only interested in leaves
                    if children:
                        continue
                    # only use categories with at least one given tag if specified
                    if len(self.category_tags) > 0 and not category.has_tag(self.category_tags, mode=any):
                        if len(hist_dict) > 0:
                            raise Exception("category {} has no required tag, but other "
                                "child categories of {} do.".format(leaf_cat, category))
                        continue

                    flavor = leaf_cat.get_aux("flavor")
                    channel = leaf_cat.get_aux("channel")
                    category_dir = input_file.GetDirectory(leaf_cat.name)

                    i_probe_jet = leaf_cat.get_aux("i_probe_jet")
                    btag_variable = btagger_cfg["variable"]

                    for process_key in category_dir.GetListOfKeys():
                        process = self.config_inst.get_process(process_key.GetName())
                        process_dir = category_dir.GetDirectory(process.name)

                        # get variable for b-tagging discriminant of probe jet
                        if process.is_data:
                            hist_shift = "nominal"
                        else:  # TODO: make nicer
                            hist_shift = self.effective_shift if not (self.iteration == 0 and not self.effective_shift.startswith("jes")) else "nominal"
                        variable_name = "jet{}_{}_{}_{}".format(i_probe_jet, btag_variable, region, hist_shift)

                        # we cannot distinguish flavors in data
                        if process.is_data and flavor != "inclusive":
                            continue
                        elif process.is_mc and flavor == "inclusive":
                            continue

                        # determine the component (heavy, light, or data) to which the flavor
                        # belongs in that region
                        if process.is_data:
                            component = "data"
                        else:
                            component = self.get_flavor_component(flavor, region)

                        # create a new hist that merges variables from multiple categories, or add
                        # to the existing one
                        hist = process_dir.Get(variable_name)

                        # rebin
                        btag_edges = self.config_inst.get_aux("binning")[region][self.b_tagger]["measurement"]
                        if self.optimize_binning:
                            binning_category = leaf_cat.get_aux("binning_category", leaf_cat)
                            btag_edges = category_binnings.get(binning_category.name,
                                btag_edges)
                        btag_edges = array.array("d", btag_edges)

                        n_bins = len(btag_edges) - 1
                        hist_rebinned = hist.Rebin(n_bins, "rebinned_{}".format(category.name), btag_edges)

                        # scale overall mc rate (per channel)
                        if process.is_mc:
                            hist_rebinned.Scale(scales[channel.name][region])

                        if component in hist_dict[category]:
                            hist_dict[category][component].Add(hist_rebinned)
                        else:
                            name = "scale_factor_{}".format(category.name)
                            hist_dict[category][component] = hist_rebinned.Clone(name)

                # calculate scale factors
                data_hist = hist_dict[category]["data"]
                lf_hist = hist_dict[category]["light"]
                hf_hist = hist_dict[category]["heavy"]

                # for the sfs, it's convenient to start with the data hist
                sf_hist = data_hist.Clone("sf_{}".format(category.name))

                # systematic shift for purity uncertainty (contamination)
                if self.effective_shift in ["lf_up", "lf_down", "hf_up", "hf_down"]:
                    contamination_scales = self.config_inst.get_aux("contamination_factors")
                    contamination_factor = contamination_scales[self.effective_shift]
                    # scale light flavour contamination in heavy flavour region
                    if self.effective_shift.split("_")[0] == "lf" and region == "hf":
                        lf_hist.Scale(contamination_factor)
                    # scale heavy flavour contamination in light flavour region
                    elif self.effective_shift.split("_")[0] == "hf" and region == "lf":
                        hf_hist.Scale(contamination_factor)

                # normalize MC histograms
                norm_factor = hist_integral(data_hist) / (hist_integral(lf_hist) + hist_integral(hf_hist))
                lf_hist.Scale(norm_factor)
                hf_hist.Scale(norm_factor)

                # subtract lf contamination from hf and vice versa
                # and do the actual division to compute scale factors
                # (this is where the physics happens)
                if region == "hf":
                    sf_hist.Add(lf_hist, -1.)
                    sf_hist.Divide(hf_hist)
                elif region == "lf":
                    sf_hist.Add(hf_hist, -1.)
                    sf_hist.Divide(lf_hist)

                # systematic shift for statistical uncertainties
                stat_uncertainties = ["{}_stats{}_{}".format(*tpl) for tpl in itertools.product(
                    ["lf", "hf"], ["1", "2"], ["up", "down"]
                )]
                if self.effective_shift in stat_uncertainties:
                    shift_flavor, shift_type, shift_direction = self.effective_shift.split("_")
                    if shift_flavor == region:
                        nbins = sf_hist.GetNbinsX()
                        for bin_idx in range(1, nbins + 1):
                            bin_center = sf_hist.GetBinCenter(bin_idx)
                            bin_content = sf_hist.GetBinContent(bin_idx)
                            bin_error = sf_hist.GetBinError(bin_idx)
                            # handle csv values smaller than 0
                            if bin_center < 0.:
                                bin_center = 0.

                            if shift_type == "stats1":
                                shift_value = bin_error * (1. - 2 * bin_center)
                            elif shift_type == "stats2":
                                shift_value = bin_error * (1. - 6 * bin_center * (1. - bin_center))
                            else:
                                raise ValueError("Unknown shift type {}".format(shift_type))
                            shift_sign = {"up": 1., "down": -1.}[shift_direction]
                            sf_hist.SetBinContent(bin_idx, bin_content + shift_sign * shift_value)

                # store the corrected sf hist
                sf_dict[category] = sf_hist

            # open the output file
            with outp["scale_factors"].localize("w") as tmp:
                with tmp.dump("RECREATE") as output_file:
                    for category in categories:
                        category_dir = output_file.mkdir(category.name)
                        category_dir.cd()
                        sf_dict[category].Write("sf")

            # for the first iteration, also save the channel rate scale factors
            if self.iteration == 0 and self.shift == "nominal":
                outp["channel_scales"].dump(scales, indent=4)


class MeasureCScaleFactors(MeasureScaleFactors):
    shifts = format_shifts(["c_stats1", "c_stats2"])

    def requires(self):
        skip_shifts = format_shifts(["hf", "lf_stats1", "lf_stats2"])
        skip_shifts = skip_shifts | jes_total_shifts

        reqs = {}
        reqs["scale_factors"] = {
            shift: MeasureScaleFactors.req(self, shift=shift,
                version=self.get_version(MeasureScaleFactors), _prefer_cli=["version"])
            for shift in MeasureScaleFactors.shifts if shift not in skip_shifts
        }
        reqs["norm"] = MergeScaleFactorWeights.req(self, normalize_cerrs=False,
                version=self.get_version(MergeScaleFactorWeights), _prefer_cli=["version"])
        if self.optimize_binning:
            reqs["binning"] = OptimizeBinning.req(self, version=self.get_version(OptimizeBinning),
                _prefer_cli=["version"])
        return reqs

    def output(self):
        return {"scale_factors": self.wlcg_target("scale_factors.root")}

    def run(self):
        import ROOT

        inp = self.input()
        outp = self.output()

        # get categories in which we measure the scale factors
        # these are stored in the config itself as we measure them inclusively over channels
        categories = []
        for category, _, _ in self.config_inst.walk_categories():
            if category.has_tag(("c", self.b_tagger), mode=all) and category.get_aux("phase_space") == "measure":
                if len(self.category_tags) > 0 and not category.has_tag(self.category_tags, mode=any):
                    continue
                categories.append(category)

        # get category-dependent binning if optimized binning is used
        if self.optimize_binning:
            category_binnings = inp["binning"].load()

        binning = self.config_inst.get_aux("binning")
        btag_edges = binning["hf"][self.b_tagger]["measurement"]
        n_bins = len(btag_edges) - 1

        # create histogram for c flavour nominal, and up and down shifts
        sf_dict = {}
        for category in categories:
            if self.optimize_binning:
                binning_category = category.get_aux("binning_category", category)
                cat_btag_edges = category_binnings.get(binning_category.name,
                    btag_edges)
            else:
                cat_btag_edges = btag_edges
            cat_btag_edges = array.array("d", cat_btag_edges)

            nominal_hist = ROOT.TH1F("sf {}".format(category.name), "Scale factors c",
                len(cat_btag_edges) - 1, cat_btag_edges
            )
            for bin_idx in range(1, len(cat_btag_edges)):
                nominal_hist.SetBinContent(bin_idx, 1.)
            sf_dict[category] = nominal_hist

        # category -> bin -> value
        sum_sqerrors_up = defaultdict(lambda: defaultdict(float))
        sum_sqerrors_down = defaultdict(lambda: defaultdict(float))

        # get nominal b flavour histograms
        nominal_hists = {}
        inp_files = inp["scale_factors"].pop("nominal")
        with inp_files["scale_factors"].load(formatter="root") as inp_file:
            # get normalization info
            norm_factors = inp["norm"].load()["nominal"]

            for category in categories:
                # get c stat errors from hf scale factors
                inp_cat_name = category.name.replace("_c_", "_hf_")

                inp_hist = inp_file.Get(inp_cat_name).Get("sf")
                inp_hist.Scale(norm_factors[inp_cat_name])

                # Decouple from file
                inp_hist.SetDirectory(0)
                nominal_hists[inp_cat_name] = inp_hist

        for shift, inp_files in inp["scale_factors"].items():
            with inp_files["scale_factors"].load(formatter="root") as inp_file:
                # get normalization info
                norm_factors = inp["norm"].load()[shift]

                for category in categories:
                    # get c stat errors from hf scale factors
                    inp_cat_name = category.name.replace("_c_", "_hf_")

                    inp_hist = inp_file.Get(inp_cat_name).Get("sf")
                    inp_hist.Scale(norm_factors[inp_cat_name])

                    # divide by nominal hist to get relative errors
                    inp_hist.Divide(nominal_hists[inp_cat_name])

                    for bin_idx in range(1, n_bins + 1):
                        error = inp_hist.GetBinContent(bin_idx) - 1.
                        if error > 0:
                            sum_sqerrors_up[category][bin_idx] += error**2
                        else:
                            sum_sqerrors_down[category][bin_idx] += error**2

        for category, sf_hist in sf_dict.items():
            # c scale factors are determined in hf categories
            inp_cat_name = category.name.replace("_c_", "_hf_")

            _, shift_type, shift_direction = self.effective_shift.split("_")
            for bin_idx in range(1, n_bins + 1):
                bin_center = sf_hist.GetBinCenter(bin_idx)

                # double the total error on the b flavour scale factors
                bin_error_up = 2 * sum_sqerrors_up[category][bin_idx]**0.5
                bin_error_down = 2 * sum_sqerrors_down[category][bin_idx]**0.5

                # compare with b scale factor, use larger difference
                inp_cat_name = category.name.replace("_c_", "_hf_")
                b_sf = nominal_hists[inp_cat_name].GetBinContent(bin_idx)
                if (b_sf - 1. > bin_error_up):
                    bin_error_up = b_sf - 1
                if (b_sf - 1. < -bin_error_down):
                    bin_error_down = abs(b_sf - 1.)
                bin_error = max([bin_error_up, bin_error_down])

                if shift_type == "stats1":
                    if bin_center < 0:
                        shift_value = bin_error_up if shift_direction == "up" else bin_error_down
                    else:
                        shift_value = bin_error * (1. - 2 * bin_center)
                elif shift_type == "stats2":
                    if bin_center < 0:
                        shift_value = 0.
                    else:
                        shift_value = bin_error * (1. - 6 * bin_center * (1. - bin_center))
                else:
                    raise ValueError("Unknown shift type {}".format(shift_type))

                shift_sign = {"up": 1., "down": -1.}[shift_direction]
                content = 1. + shift_sign * shift_value
                sf_hist.SetBinContent(bin_idx, content if content > 0. else 0.)

        # open the output file
        with outp["scale_factors"].localize("w") as tmp:
            with tmp.dump("RECREATE") as output_file:
                for category in categories:
                    category_dir = output_file.mkdir(category.name)
                    category_dir.cd()
                    sf_dict[category].Write("sf")


class MeasureCScaleFactorsWrapper(WrapperTask):

    wrapped_task = MeasureCScaleFactors


class FitScaleFactors(MeasureScaleFactors):

    fix_normalization = luigi.BoolParameter()

    shifts = MeasureScaleFactors.shifts | MeasureCScaleFactors.shifts

    def __init__(self, *args, **kwargs):
        super(FitScaleFactors, self).__init__(*args, **kwargs)

        self.has_c_shift = self.shift in MeasureCScaleFactors.shifts

    def requires(self):
        reqs = {}
        # get scale factor histograms
        if self.has_c_shift:
            reqs["sf"] = MeasureCScaleFactors.req(self,
                version=self.get_version(MeasureCScaleFactors), _prefer_cli=["version"])
        else:
            reqs["sf"] = MeasureScaleFactors.req(self,
                version=self.get_version(MeasureScaleFactors), _prefer_cli=["version"])

        # get scaling factors to normalize scale factors
        if self.fix_normalization:
            if self.has_c_shift:
                reqs["norm"] = MergeScaleFactorWeights.req(self, normalize_cerrs=True,
                    version=self.get_version(MergeScaleFactorWeights), _prefer_cli=["version"])
            else:
                reqs["norm"] = MergeScaleFactorWeights.req(self, normalize_cerrs=False,
                version=self.get_version(MergeScaleFactorWeights), _prefer_cli=["version"])
        return reqs

    def store_parts(self):
        normalization_part = "rescaled" if self.fix_normalization else "unscaled"
        return super(FitScaleFactors, self).store_parts() + (normalization_part,)

    def output(self):
        outp = {
            "sf": self.wlcg_target("scale_factors.root")
        }
        # if the scale factors are normalized, this is the final task
        # and the .csv output for the b-tag reader should be created
        if self.fix_normalization:
            outp["csv"] = self.wlcg_target("scale_factors.csv")
            outp["functions"] = self.wlcg_target("sf_functions.root")
        return outp

    def run(self):
        import ROOT

        inp = self.input()
        outp = self.output()

        interpolation_bins = 1000

        # cannot get the function from ROOT, use scipy instead
        from scipy.interpolate import PchipInterpolator, BPoly, PPoly

        # get categories in which to fit the scale factors
        categories = []
        for category, _, _ in self.config_inst.walk_categories():
            if len(self.category_tags) > 0 and not category.has_tag(self.category_tags, mode=any):
                continue
            if self.has_c_shift:
                if category.get_aux("region", None) == "c" and category.get_aux("phase_space") == "measure":
                    if category.has_tag(self.b_tagger):
                        categories.append(category)
            else:
                if category.has_tag(("merged", self.b_tagger), mode=all) and category.get_aux("phase_space") == "measure":
                    categories.append(category)

        # get scaling factors for normalization
        if self.fix_normalization:
            norm_factors = inp["norm"].load()[self.effective_shift]

        # contents of .csv file for scale factors
        fit_results = []
        # finely binned histograms to write to the output file
        hist_dict = {}
        # TF1's to write to the output file
        function_dict = {}
        with inp["sf"]["scale_factors"].load("r") as input_file:
            for category in categories:
                region = category.get_aux("region")
                category_dir = input_file.GetDirectory(category.name)

                # get scale factor histogram
                hist_keys = category_dir.GetListOfKeys()
                if len(hist_keys) != 1:
                    raise ValueError("Found more than one histogram in %s, cannot identify scale "
                        "factor hist." % category_dir)
                hist = category_dir.Get(hist_keys[0].GetName())
                nbins = hist.GetNbinsX()
                if self.fix_normalization:
                    hist.Scale(norm_factors[category.name])

                x_axis = hist.GetXaxis()
                interpolation_hist = ROOT.TH1D(hist.GetName() + "_fine", hist.GetTitle(),
                    interpolation_bins, x_axis.GetXmin(), x_axis.GetXmax())
                x_values = ROOT.vector("double")()
                y_values = ROOT.vector("double")()
                for bin_idx in range(1, nbins + 1):
                    if hist.GetBinCenter(bin_idx) < 0:
                        continue
                    x_values.push_back(hist.GetBinCenter(bin_idx))
                    y_values.push_back(hist.GetBinContent(bin_idx))

                interpolator = PchipInterpolator(x_values, y_values)
                # define region in which to use interpolation
                first_point, last_point = min(x_values), max(x_values)

                # create finely binned histogram from either TF1 or interpolator
                for bin_idx in range(interpolation_bins + 2):
                    bin_center = interpolation_hist.GetBinCenter(bin_idx)
                    if bin_center < 0:
                        interpolation_hist.SetBinContent(bin_idx, hist.GetBinContent(1))
                    elif bin_center < first_point:
                        interpolation_hist.SetBinContent(bin_idx, interpolator(first_point))
                    elif bin_center > last_point:
                        interpolation_hist.SetBinContent(bin_idx, interpolator(last_point))
                    else:
                        interpolation_hist.SetBinContent(bin_idx, interpolator(bin_center))
                hist_dict[category] = interpolation_hist

                # fill .csv file in final iteration (after normalization fix)
                # also create piecewise linear TF1's
                if self.fix_normalization:
                    function_pieces = []

                    results = {}
                    results["eta_min"], results["eta_max"] = category.get_aux("eta")
                    pt_range = category.get_aux("pt")
                    results["pt_min"] = pt_range[0]
                    results["pt_max"] = min(pt_range[1], 10000.)  # replace inf
                    results["flavor_id"] = self.config_inst.get_aux("flavor_ids")[region]

                    if self.effective_shift == "nominal":
                        sysType = "central"
                    else:
                        sys_name, direction = self.effective_shift.rsplit("_", 1)
                        sysType = "{}_{}".format(direction,
                            sys_name.replace("c_stats", "cferr").replace("lf_stats", "lfstats").replace("hf_stats", "hfstats")
                        )
                    results["sysType"] = sysType

                    # skip unwanted combinations
                    if "cferr" in sysType and region != "c":
                        continue

                    fit_results_tpl = "3, iterativefit, {sysType}, {flavor_id}, {eta_min}, " \
                        "{eta_max}, {pt_min}, {pt_max}".format(**results)
                    fit_results.append(fit_results_tpl + ", -15, 0, {}".format(hist.GetBinContent(1)))
                    fit_results.append(fit_results_tpl + ", 0, {}, {}".format(first_point,
                        interpolator(first_point)))

                    function_pieces.append("(x < 0) * {}".format(hist.GetBinContent(1)))
                    function_pieces.append("(x >= 0) * (x < {}) * {}".format(first_point,
                        interpolator(first_point)))

                    # intermediate functions
                    # change interpolated function from bernstein to power basis
                    bpoly_interpolation = BPoly(interpolator.c, interpolator.x)
                    ppoly_interpolation = PPoly.from_bernstein_basis(bpoly_interpolation)

                    interpolator_idx = 0
                    for bin_idx in range(1, nbins):
                        if hist.GetBinCenter(bin_idx) < first_point:
                            continue
                        x_min = hist.GetBinCenter(bin_idx)
                        x_max = hist.GetBinCenter(bin_idx + 1)

                        interpolator_coefficients = ppoly_interpolation.c[:, interpolator_idx]
                        interpolator_x = ppoly_interpolation.x[interpolator_idx]
                        interpolator_idx += 1
                        formula = ""
                        for i in xrange(3):
                            formula += "{}*".format(interpolator_coefficients[i]) + \
                                "*".join(["(x-{})".format(interpolator_x)] * (3 - i))
                            formula += "+" if interpolator_coefficients[i + 1] >= 0. else ""
                        formula += str(interpolator_coefficients[3])
                        fit_results.append(fit_results_tpl + ", {}, {}, {}".format(x_min,
                            x_max, formula))
                        function_pieces.append("(x >= {}) * (x < {}) * ({})".format(
                            x_min, x_max, formula))

                    fit_results.append(fit_results_tpl + ", {}, 1.1, {}".format(last_point,
                        interpolator(last_point)))
                    function_pieces.append("(x >= {}) * {}".format(last_point,
                        interpolator(last_point)))

                    function = ROOT.TF1("sf_{}".format(category.name),
                        " + ".join(["({})".format(piece) for piece in function_pieces]),
                        -2., 1.1)
                    function_dict[category] = function

                    # sanity check
                    for val in x_values:
                        if (abs(function.Eval(val) - interpolator(val)) > 1e-5):
                            raise Exception("SF Function does not match interpolator values: "
                                "{} vs {}".format(function.Eval(val), interpolator(val)))

            # write to output file
            with outp["sf"].localize("w") as tmp:
                with tmp.dump("RECREATE") as output_file:
                    for category, hist in hist_dict.items():
                        category_dir = output_file.mkdir(category.name)
                        category_dir.cd()
                        hist.Write("sf")
            if self.fix_normalization:
                with outp["csv"].localize("w") as tmp:
                    with tmp.open("w") as result_file:
                        result_file.write("\n".join(fit_results))

                with outp["functions"].localize("w") as tmp:
                    with tmp.dump("RECREATE") as output_file:
                        for category, func in function_dict.items():
                            category_dir = output_file.mkdir(category.name)
                            category_dir.cd()
                            func.Write("sf")


class FitScaleFactorsWrapper(WrapperTask):

    wrapped_task = FitScaleFactors


# bundle scale factors to reduce number of dCache requests
class BundleScaleFactors(AnalysisTask):

    iteration = FitScaleFactors.iteration
    b_tagger = FitScaleFactors.b_tagger
    optimize_binning = FitScaleFactors.optimize_binning
    category_tags = FitScaleFactors.category_tags

    fix_normalization = FitScaleFactors.fix_normalization
    include_cshifts = luigi.BoolParameter()

    shifts = MeasureScaleFactors.shifts | MeasureCScaleFactors.shifts

    def __init__(self, *args, **kwargs):
        super(BundleScaleFactors, self).__init__(*args, **kwargs)

        self.shifts = MeasureScaleFactors.shifts
        if self.include_cshifts:
            self.shifts = self.shifts | MeasureCScaleFactors.shifts

    def requires(self):
        reqs = {
            shift: FitScaleFactors.req(self, shift=shift,
            version=self.get_version(FitScaleFactors), _prefer_cli=["version"])
            for shift in self.shifts
        }
        return reqs

    def store_parts(self):
        shift_part = "all" if self.include_cshifts else "no_cshift"
        normalization_part = "rescaled" if self.fix_normalization else "unscaled"
        binning_part = "optimized" if self.optimize_binning else "default"
        return super(BundleScaleFactors, self).store_parts() + (self.b_tagger,) \
            + (self.iteration,) + (normalization_part,) + (shift_part,) + (binning_part,)

    def output(self):
        return self.wlcg_target("scale_factors.root")

    def run(self):
        import ROOT

        inp = self.input()
        outp = self.output()

        with outp.localize("w") as tmp_out:
            with tmp_out.dump("RECREATE") as output_file:
                for shift, input_targets in inp.items():
                    shift_dir = output_file.mkdir(shift)

                    with input_targets["sf"].load("r") as input_file:
                        for category in input_file.GetListOfKeys():
                            category_dir = shift_dir.mkdir(category.GetName())
                            category_dir.cd()
                            hist = input_file.Get(category.GetName()).Get("sf")
                            hist.Write("sf")


class CreateScaleFactorResults(AnalysisTask):

    iteration = MeasureScaleFactors.iteration
    b_tagger = MeasureScaleFactors.b_tagger
    optimize_binning = MeasureScaleFactors.optimize_binning

    def requires(self):
        reqs = {shift: FitScaleFactors.req(self, shift=shift, fix_normalization=True,
            version=self.get_version(FitScaleFactors), _prefer_cli=["version"])
            for shift in FitScaleFactors.shifts}
        return reqs

    def store_parts(self):
        binning_part = "optimized" if self.optimize_binning else "default"
        return super(CreateScaleFactorResults, self).store_parts() + (self.b_tagger,) + (self.iteration,) \
            + (binning_part,)

    def output(self):
        outp = {}
        outp["csv"] = self.wlcg_target("scale_factors.csv")
        outp["root_lf"] = self.wlcg_target("scale_factors_lf.root")
        outp["root_hf"] = self.wlcg_target("scale_factors_hf.root")
        return outp

    def run(self):
        def get_func_name(category, shift):
            _, region, pt_range, eta_range, _ = category.split("__")
            name = "c_csv" if region == "c" else "csv"
            name += "_ratio"

            # get pt and eta indices of category
            pt_lower = int(pt_range.split("To")[0][2:])
            eta_lower = eta_range.split("To")[0][3:]

            if region in ("hf", "c"):
                eta_bin = 0
                pt_bin = {20: 0, 30: 1, 50: 2, 70: 3, 100: 4, 140: 5}[pt_lower]
            elif region == "lf":
                eta_bin = {"0": 0, "0p8": 1, "1p6": 2}[eta_lower]
                pt_bin = {20: 0, 30: 1, 40: 2, 60: 3, 100: 4}[pt_lower]
            else:
                raise ValueError("Unknown region {}".format(region))

            name += "_Pt{}".format(pt_bin)
            name += "_Eta{}".format(eta_bin)

            name += "_final"

            # Parse uncertainty names
            if shift != "nominal":
                uncertainty = shift.replace("_", "").replace("lf", "LF").replace("hf", "HF")
                uncertainty = uncertainty.replace("jes", "JES").replace("up", "Up").replace("down", "Down")

                if "cstats" in uncertainty and region == "c":
                    uncertainty = uncertainty.replace("cstats", "cErr")
                elif "cstats" in uncertainty:
                    return None, None

                name += "_{}".format(uncertainty)
            return name, region

        import ROOT

        inp = self.input()
        outp = self.output()

        csv_results = ["3, iterativefit, central, 1, {}, {}, 20.0, 10000, -15, 1.1, 1.0\n".format(*self.config_inst.get_aux("binning")["c"]["abs(eta)"])]
        hf_funcs = []
        lf_funcs = []

        for shift, inp_files in inp.items():
            # combine csv files for b-tag reader
            with inp_files["csv"].open("r") as csv_file:
                for line in csv_file.readlines():
                    if not line.endswith("\n"):
                        line += "\n"
                    line = line.replace("jesTotal", "jes")
                    csv_results.append(line)

            # combine functions files
            with inp_files["functions"].load("r") as root_file:
                for category_key in root_file.GetListOfKeys():
                    category_dir = root_file.Get(category_key.GetName())
                    func = category_dir.Get("sf")
                    name, region = get_func_name(category_key.GetName(), shift)
                    if name is None:
                        continue

                    if region in ("hf", "c"):
                        hf_funcs.append((name, func))
                    elif region == "lf":
                        lf_funcs.append((name, func))
                    else:
                        raise ValueError("Unknown region {}".format(region))

        # add nominal c tag function (flat value of 1)
        c_func = ROOT.TF1("c nominal", "1.", -2., 1.1)
        for iPt in range(6):
            hf_funcs.append(("c_csv_ratio_Pt{}_Eta0_final".format(iPt), c_func.Clone()))

        with outp["root_lf"].localize("w") as tmp:
            with tmp.dump("RECREATE") as output_file:
                for name, func in lf_funcs:
                    func.Write(name)

        with outp["root_hf"].localize("w") as tmp:
            with tmp.dump("RECREATE") as output_file:
                for name, func in hf_funcs:
                    func.Write(name)

        with outp["csv"].localize("w") as tmp:
            with tmp.open("w") as result_file:
                result_file.write("".join(csv_results))
