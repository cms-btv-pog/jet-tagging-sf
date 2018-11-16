# -*- coding: utf-8 -*-


import os
import array
import numpy as np

from collections import itertools, defaultdict

from analysis.tasks.base import ShiftTask
from analysis.tasks.hists import MergeHistograms


class MeasureScaleFactors(ShiftTask):

    iteration = MergeHistograms.iteration

    shifts = MergeHistograms.shifts

    def requires(self):
        reqs = {
            "hist": MergeHistograms.req(self, branch=0, version=self.get_version(MergeHistograms),
                _prefer_cli=["version"])
        }
        if self.iteration > 0:
            reqs["scale"] = MeasureScaleFactors.req(self, iteration=0,
                version=self.get_version(MeasureScaleFactors), _prefer_cli=["version"])
        return reqs

    def store_parts(self):
        return super(MeasureScaleFactors, self).store_parts() + (self.iteration,)

    def output(self):
        outputs = {"scale_factors": self.wlcg_target("scale_factors.root")}
        if self.iteration == 0:
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

        import ROOT

        inp = self.input()
        outp = self.output()

        # get categories in which we measure the scale factors
        # these are stored in the config itself as we measure them inclusively over channels
        categories = []
        for category, _, _ in self.config_inst.walk_categories():
            if category.has_tag("merged") and category.get_aux("phase_space") == "measure":
                categories.append(category)

        # get categories from which to determine the rate scaling of MC to data
        # only needed for the first iteration, where the scaling is saved for further use
        if self.iteration == 0:
            scale_categories = {}
            for channel in self.config_inst.channels:
                scale_categories[channel] = {}
                for category, _, children in channel.walk_categories():
                    if category.has_tag("scales") and category.get_aux("phase_space") == "closure":
                        region = category.get_aux("region")
                        scale_categories[channel][region] = category

        btagger_cfg = self.config_inst.get_aux("btagger")

        # category -> component (heavy/light) -> histogram
        hist_dict = {}
        # category -> histogram
        sf_dict = {}

        # create n-d histograms to hold all scale factors for lf/hf jets
        binning = self.config_inst.get_aux("binning")
        sf_hists_nd = {}
        for region in ["lf", "hf"]:
            eta_edges = array.array("d", binning[region]["abs(eta)"])
            pt_edges = array.array("d", binning[region]["pt"])
            btag_edges = array.array("d", binning[region][btagger_cfg["name"]]["measurement"])
            sf_hist = ROOT.TH3F(
                "scale_factors_{}".format(region), "Scale factors {}".format(region),
                len(eta_edges) - 1, eta_edges,
                len(pt_edges) - 1, pt_edges,
                len(btag_edges) - 1, btag_edges,
            )
            sf_hists_nd[region] = sf_hist

        with inp["hist"].load("r") as input_file:
            # get scale factor to scale MC (withouts b-tag SFs) to data per channel
            if self.iteration == 0:
                variable_name = "jet1_pt"
                scales = defaultdict(dict)
                for channel, region_categories in scale_categories.items():
                    for region, category in region_categories.items():
                        data_yield, mc_yield = 0., 0.
                        category_dir = input_file.GetDirectory(category.name)

                        # sum over processes to get mc and data yields
                        for process_key in category_dir.GetListOfKeys():
                            process = self.config_inst.get_process(process_key.GetName())
                            process_dir = category_dir.GetDirectory(process.name)

                            hist = process_dir.Get(variable_name)
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

                    flavor = leaf_cat.get_aux("flavor")
                    channel = leaf_cat.get_aux("channel")
                    category_dir = input_file.GetDirectory(leaf_cat.name)

                    # get variable for b-tagging discriminant of probe jet
                    i_probe_jet = leaf_cat.get_aux("i_probe_jet")
                    btag_variable = btagger_cfg["variable"]
                    variable_name = "jet{}_{}_{}".format(i_probe_jet, btag_variable, region)

                    for process_key in category_dir.GetListOfKeys():
                        process = self.config_inst.get_process(process_key.GetName())
                        process_dir = category_dir.GetDirectory(process.name)

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
                        btag_edges = array.array("d", binning[region][btagger_cfg["name"]]["measurement"])
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
                if self.shift in ["lf_up", "lf_down", "hf_up", "hf_down"]:
                    contamination_scales = self.config_inst.get_aux("contamination_factors")
                    contamination_factor = contamination_scales[self.shift]
                    # scale light flavour contamination in heavy flavour region
                    if self.shift.split("_")[0] == "lf" and region == "hf":
                        lf_hist.Scale(contamination_factor)
                    # scale heavy flavour contamination in light flavour region
                elif self.shift.split("_")[0] == "hf" and region == "lf":
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
                if self.shift in stat_uncertainties:
                    shift_flavor, shift_type, shift_direction = self.shift.split("_")
                    if shift_flavor == region:
                        nbins = sf_hist.GetNbinsX()
                        for bin_idx in range(1, nbins + 1):
                            bin_center = sf_hist.GetBinCenter(bin_idx)
                            bin_content = sf_hist.GetBinContent(bin_idx)
                            bin_error = sf_hist.GetBinError(bin_idx)
                            if shift_type == "stats1":
                                shift_value = bin_error * (1. - 2 * bin_center)# TODO: Bin at CSV < 0
                            elif shift_type == "stats2":
                                shift_value = bin_error * (1. - 6 * bin_center * (1. - bin_center))
                            else:
                                raise ValueError("Unknown shift type {}".format(shift_type))
                            shift_sign = {"up": 1., "down": -1.}[shift_direction]
                            sf_hist.SetBinContent(bin_content + shift_sign * shift_value)

                # store the corrected sf hist
                sf_dict[category] = sf_hist

                # write to nd histograms
                eta_val = np.mean(category.get_aux("eta"))
                pt_val = np.mean(category.get_aux("pt"))
                for bin_idx in range(1, sf_hist.GetNbinsX() + 1):  # TODO: Add over- and underflow bin
                    sf_hists_nd[region].Fill(eta_val, pt_val, sf_hist.GetBinCenter(bin_idx),
                        sf_hist.GetBinContent(bin_idx))

            # open the output file
            with outp["scale_factors"].localize("w") as tmp:
                with tmp.dump("RECREATE") as output_file:
                    for category in categories:
                        category_dir = output_file.mkdir(category.name)
                        category_dir.cd()
                        sf_dict[category].Write("sf")
                    for region in sf_hists_nd:
                        region_dir = output_file.mkdir(region)
                        region_dir.cd()
                        sf_hists_nd[region].Write("sf")

            # for the first iteration, also save the channel rate scale factors
            if self.iteration == 0:
                outp["channel_scales"].dump(scales, indent=4)


class FitScaleFactors(MeasureScaleFactors):

    final_it = MergeHistograms.final_it

    def requires(self):
        return MeasureScaleFactors.req(self, version=self.get_version(MeasureScaleFactors),
            _prefer_cli=["version"])

    def output(self):
        outp = {
            "sf": self.wlcg_target("scale_factors.root")
        }
        if self.final_it:
            outp["csv"] = self.wlcg_target("scale_factors.csv")
        return outp

    def run(self):
        import ROOT

        inp = self.input()
        outp = self.output()

        interpolation_type = ROOT.Math.Interpolation.kAKIMA
        interpolation_bins = 1000

        def fit_func_pol6(x_min=0.0, x_max=0.5):  # 6th degree polynomial for LF region
            func = ROOT.TF1("f_pol6", "pol6(0)", x_min, x_max)
            for i_param in range(func.GetNpar()):
                func.SetParameter(i_param, 1.)
            return func

        def fit_func_pol1(x_min, x_max, init_params=[1., 1.]):  # 1st degree polynomial for final fit in HF region
            func = ROOT.TF1("f_pol1", "pol1(0)", x_min, x_max)
            if not len(init_params) == func.GetNpar():
                raise ValueError("Expected {} parameter values, but got {}".format(func.GetNpar(),
                    init_params))
            func.SetParameters(*init_params)
            return func

        # get categories in which to fit the scale factors
        categories = []
        for category, _, _ in self.config_inst.walk_categories():
            if category.has_tag("merged") and category.get_aux("phase_space") == "measure":
                categories.append(category)

        # contents of .csv file for scale factors
        fit_results = []
        # finely binned histograms to write to the output file
        hist_dict = {}
        with inp["scale_factors"].load("r") as input_file:
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

                x_axis = hist.GetXaxis()
                interpolation_hist = ROOT.TH1D(hist.GetName() + "_fine", hist.GetTitle(),
                    interpolation_bins, x_axis.GetXmin(), x_axis.GetXmax())
                if region == "lf":
                    fit_function = fit_func_pol6()
                    # perform fit
                    hist.Fit(fit_function, "+mrNQ0S")
                    interpolator = fit_function
                    # define region in which to use the function values to create the histogram
                    # (centers of second and second to last bin)
                    first_point = hist.GetBinCenter(2)
                    last_point = hist.GetBinCenter(nbins - 1)
                elif region == "hf":
                    x_values = ROOT.vector("double")()
                    y_values = ROOT.vector("double")()
                    for bin_idx in range(1, nbins + 1):
                        if hist.GetBinCenter(bin_idx) < 0:
                            continue
                        x_values.push_back(hist.GetBinCenter(bin_idx))
                        y_values.push_back(hist.GetBinContent(bin_idx))
                    interpolator = ROOT.Math.Interpolator(x_values, y_values, interpolation_type)
                    # define region in which to use interpolation
                    first_point, last_point = min(x_values), max(x_values)
                else:
                    raise ValueError("Unknown region %s" % region)

                # create finely binned histogram from either TF1 or interpolator
                for bin_idx in range(interpolation_bins + 2):
                    bin_center = interpolation_hist.GetBinCenter(bin_idx)
                    if bin_center < 0:
                        interpolation_hist.SetBinContent(bin_idx, hist.GetBinContent(1))
                    elif bin_center < first_point:
                        interpolation_hist.SetBinContent(bin_idx, interpolator.Eval(first_point))
                    elif bin_center > last_point:
                        interpolation_hist.SetBinContent(bin_idx, interpolator.Eval(last_point))
                    else:
                        interpolation_hist.SetBinContent(bin_idx, interpolator.Eval(bin_center))
                hist_dict[category] = interpolation_hist

                # fill .csv file in final iteration
                if self.final_it:
                    results = {}
                    results["eta_min"], results["eta_max"] = category.get_aux("eta")
                    pt_range = category.get_aux("pt")
                    results["pt_min"] = pt_range[0]
                    results["pt_max"] = min(pt_range[1], 10000.)  # replace inf
                    results["flavor_id"] = self.config_inst.get_aux("flavor_ids")[region]

                    fit_results_tpl = "3, iterativefit, central, {flavor_id}, {eta_min}, " \
                        "{eta_max}, {pt_min}, {pt_max}".format(**results)
                    fit_results.append(fit_results_tpl + ", -15, 0, {}".format(hist.GetBinContent(1)))
                    fit_results.append(fit_results_tpl + ", 0, {}, {}".format(first_point,
                        interpolator.Eval(first_point)))

                    # intermediate functions
                    if region == "lf":
                        fit_results.append(fit_results_tpl + ", {}, {}, {}".format(first_point,
                            last_point, str(interpolator.GetExpFormula("p"))))
                    elif region == "hf":  # piecewise linear function
                        for bin_idx in range(1, nbins):
                            if hist.GetBinCenter(bin_idx) < first_point:
                                continue
                            x_min = hist.GetBinCenter(bin_idx)
                            x_max = hist.GetBinCenter(bin_idx + 1)
                            y_start = hist.GetBinContent(bin_idx)
                            y_end = hist.GetBinContent(bin_idx + 1)
                            slope = (y_end - y_start) / (x_max - x_min)
                            intercept = y_start - x_min * slope
                            func = fit_func_pol1(x_min, x_max, [intercept, slope])

                            fit_results.append(fit_results_tpl + ", {}, {}, {}".format(x_min,
                                x_max, str(func.GetExpFormula("p"))))
                    else:
                        raise ValueError("Unknown region %s" % region)

                    fit_results.append(fit_results_tpl + ", {}, 1.1, {}".format(last_point,
                        interpolator.Eval(last_point)))
            # write to output file
            with outp["sf"].localize("w") as tmp:
                with tmp.dump("RECREATE") as output_file:
                    for category, hist in hist_dict.items():
                        category_dir = output_file.mkdir(category.name)
                        category_dir.cd()
                        hist.Write("sf")
            if self.final_it:
                with outp["csv"].localize("w") as tmp:
                    with tmp.open("w") as result_file:
                        result_file.write("\n".join(fit_results))
