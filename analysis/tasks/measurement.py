# -*- coding: utf-8 -*-


import os

from analysis.tasks.base import AnalysisTask
from analysis.tasks.hists import MergeHistograms


class CalculateScaleFactors(AnalysisTask):
    def requires(self):
        return MergeHistograms.req(self, branch=0, version=self.get_version(MergeHistograms),
                                   _prefer_cli=["version"])

    def output(self):
        return self.wlcg_target("sfs.root")

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

        # category -> component (heavy/light) -> histogram
        hist_dict = {}
        # category -> histogram
        sf_dict = {}

        with inp.load("r") as input_file:
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
                    btag_variable = self.config_inst.get_aux("btagger")["variable"]
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
                        if component in hist_dict[category]:
                            hist_dict[category][component].Add(hist)
                        else:
                            name = "scale_factor_{}".format(category.name)
                            hist_dict[category][component] = hist.Clone(name)

                data_hist = hist_dict[category]["data"]
                lf_hist = hist_dict[category]["light"]
                hf_hist = hist_dict[category]["heavy"]

                # for the sfs, it's convenient to start with the data hist
                sf_hist = data_hist.Clone("sf_{}".format(category.name))

                # scale overall rate of mc in this category to data
                scale = data_hist.Integral() / (lf_hist.Integral() + hf_hist.Integral())
                lf_hist.Scale(scale)
                hf_hist.Scale(scale)

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
            with outp.localize("w") as tmp:
                with tmp.dump("RECREATE") as output_file:
                    for category in categories:
                        category_dir = output_file.mkdir(category.name)
                        category_dir.cd()
                        sf_dict[category].Write("sf")
