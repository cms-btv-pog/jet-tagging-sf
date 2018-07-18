# -*- coding: utf-8 -*-


import os
import array
import numpy as np

from collections import defaultdict

from analysis.tasks.base import AnalysisTask
from analysis.tasks.hists import MergeHistograms


class MeasureScaleFactors(AnalysisTask):

    iteration = MergeHistograms.iteration

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
        if region == "HF":
            return "heavy" if flavor == "b" else "light"
        elif region == "LF":
            return "light" if flavor == "udsg" else "heavy"
        else:
            raise ValueError("unexpected region %s" % region)

    def run(self):
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
        for region in ["LF", "HF"]:
            eta_edges = array.array("f", binning[region]["eta"])
            pt_edges = array.array("f", binning[region]["pt"])
            btag_edges = array.array("f", binning[region][btagger_cfg["name"]])
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
                                data_yield += hist.Integral(0, hist.GetNbinsX() + 1)
                            else:
                                mc_yield += hist.Integral(0, hist.GetNbinsX() + 1)
                        scale = data_yield / mc_yield
                        scales[channel.name][region] = scale
            else:
                scales = inp["channel_scales"].load()
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
                        # scale overall mc rate (per channel)
                        if process.is_mc:
                            hist.Scale(scales[channel.name][region])

                        if component in hist_dict[category]:
                            hist_dict[category][component].Add(hist)
                        else:
                            name = "scale_factor_{}".format(category.name)
                            hist_dict[category][component] = hist.Clone(name)

                # calculate scale factors
                data_hist = hist_dict[category]["data"]
                lf_hist = hist_dict[category]["light"]
                hf_hist = hist_dict[category]["heavy"]

                # for the sfs, it's convenient to start with the data hist
                sf_hist = data_hist.Clone("sf_{}".format(category.name))

                # subtract lf contamination from hf and vice versa
                # and do the actual division to compute scale factors
                # (this is where the physics happens)
                if region == "HF":
                    sf_hist.Add(lf_hist, -1.)
                    sf_hist.Divide(hf_hist)
                elif region == "LF":
                    sf_hist.Add(hf_hist, -1.)
                    sf_hist.Divide(lf_hist)

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
