# -*- coding: utf-8 -*-


import os
import array

from analysis.tasks.base import AnalysisTask
from analysis.tasks.hists import MergeHistograms


class CalculateScaleFactors(AnalysisTask):

    iteration = MergeHistograms.iteration

    def requires(self):
        reqs = {
            "hist": MergeHistograms.req(self, branch=0, version=self.get_version(MergeHistograms),
                _prefer_cli=["version"])
        }
        if self.iteration > 0:
            reqs["scale"] = CalculateScaleFactors.req(self, iteration=0,
                version=self.get_version(CalculateScaleFactors), _prefer_cli=["version"])

    def store_parts(self):
        return super(CalculateScaleFactors, self).store_parts() + (self.iteration,)

    def output(self):
        outp = {
            "sfs": self.wlcg_target("sfs.root")
        }
        if self.iteration == 0:
            outp["scale"] = self.wlcg_target("scale.json")

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
            raise ValueError("Unexpected region %s" % region)

    def run(self):
        import ROOT

        inp = self.input()
        outp = self.output()
        outp.parent.touch(0o0770)

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
            for channel in self.config_inst.channels.values():
                for category, _, _ in channel.walk_categories():
                    if category.has_tag("measure"):
                        if channel in scale_categories:
                            raise Exception("Should only have one category with tag "
                            "'measure' per channel, but got "
                            "{} and {}".format(category, scale_categories[channel]))
                        scale_categories[channel] = category

        btagger_cfg = self.config_inst.get_aux("btagger")

        # category -> component (heavy/light) -> histogram
        hist_dict = {}
        # category -> histogram
        sf_dict = {}

        # create n-d histograms to hold all scale factors for lf/hf jets
        binning = self.config_inst.get_aux("binning")
        sf_hists_nd = {}
        for region in ["LF", "HF"]:
            eta_bins = array.array(binning[region]["eta"])
            pt_bins = array.array(binning[region]["eta"])
            btag_bins = array.array(binning[region][btagger_cfg["name"]])
            sf_hist = ROOT.TH3F("scale_factors_{}".format(region), "Scale factors {}".format(region),
                len(eta_bins), eta_bins, len(pt_bins), pt_bins, len(btag_bins), btag_bins)
            sf_hists_nd[region] = sf_hist

        with inp["hist"].load("r") as input_file:
            # get scale factor to scale MC (withouts b-tag SFs) to data
            if self.iteration == 0:
                variable_name = "jet1_pt"
                scales = {}
                for channel, category in scale_categories.items():
                    data_yield, mc_yield = 0., 0.
                    category_dir = input_file.GetDirectory(category.name)

                    # sum over processes to get mc and data yields
                    for process_key in category_dir.GetListOfKeys:
                        process = self.config_inst.get_process(process_key.GetName())
                        process_dir = category_dir.GetDirectory(process.name)

                        hist = process_dir.Get(variable_name)
                        if process.is_data:
                            data_yield += hist.Integral()
                        else:
                            mc_yield += hist.Integral()
                    scale = data_yield / mc_yield
                    scales[category] = scale
            else:
                scales = inp["scale"].load()

            for category in categories:
                region = category.get_aux("region")

                hist_dict[category] = {}
                for leaf_cat, _, children in category.walk_categories():
                    # we are only interested in leaves
                    if children:
                        continue

                    flavor = leaf_cat.get_aux("flavor")
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
                            hist.Scale(scales[channel.name])

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
                    sf_hist.Add(lf_hist, -1)
                    sf_hist.Divide(hf_hist)
                elif region == "LF":
                    sf_hist.Add(hf_hist, -1)
                    sf_hist.Divide(lf_hist)

                # store the corrected sf hist
                sf_dict[category] = sf_hist

            # open the output file
            with outp["sfs"].localize("w") as tmp:
                with tmp.dump("RECREATE") as output_file:
                    for category in categories:
                        category_dir = output_file.mkdir(category.name)
                        category_dir.cd()
                        sf_dict[category].Write("sf")
            # for the first iteration, save the rate scale factors
            if self.iteration == 0:
                outp["scale"].dump(scales, indent=4)
