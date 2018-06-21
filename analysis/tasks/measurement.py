# -*- coding: utf-8 -*-

import os
import re #TODO: Remove for next iteration

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
        Depending on the region, c jets are counted for either the heavy or the light flavor component.
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

        # get categories
        categories = []
        channels = self.config_inst.channels.values()

        for channel in channels:
            for category, _, _ in channel.walk_categories():
                if category.has_tag("merged"):
                    categories.append(category)

        hist_dict = {}
        sf_dict = {}
        with inp.load("r") as input_file:
            for category in categories:
                region = category.get_aux("region")

                hist_dict[category] = {}
                for child_cat, _, _ in category.walk_categories():
                    child_cat_name = re.sub(r"__f", "__", child_cat.name) # TODO: Remove for next iteration of histograms
                    category_dir = input_file.GetDirectory(child_cat_name)

                    flavor = child_cat.get_aux("flavor")

                    # get variable for b-tagging discriminant of probe jet
                    i_probe_jet = child_cat.get_aux("i_probe_jet")
                    btag_variable = self.config_inst.get_aux("btagger")["variable"]
                    variable = "jet{}_{}".format(i_probe_jet, btag_variable)

                    for process_key in category_dir.GetListOfKeys():
                        process = self.config_inst.get_process(process_key.GetName())
                        process_dir = category_dir.GetDirectory(process.name)

                        if process.is_data and flavor != "inclusive":
                            continue
                        elif process.is_mc and flavor == "inclusive":
                            continue

                        hist = process_dir.Get(variable)
                        if process.is_data:
                            component = "data"
                        else:
                            component = self.get_flavor_component(flavor, region)

                        if component in hist_dict[category]:
                            hist_dict[category][component].Add(hist)
                        else:
                            hist_dict[category][component] = hist.Clone("scale_factor_{}".format(region))

                sf_hist = hist_dict[category]["data"]
                lf_hist = hist_dict[category]["light"]
                hf_hist = hist_dict[category]["heavy"]
                # scale overall rate of mc in this category to data
                scale = sf_hist.Integral()/(lf_hist.Integral() + hf_hist.Integral())
                lf_hist.Scale(scale)
                hf_hist.Scale(scale)

                if region == "HF":
                    sf_hist.Add(lf_hist, -1)
                    sf_hist.Divide(hf_hist)
                elif region == "LF":
                    sf_hist.Add(hf_hist, -1)
                    sf_hist.Divide(lf_hist)
                sf_dict[category] = sf_hist

            # open the output file
            with outp.localize("w") as tmp:
                with tmp.dump("RECREATE") as output_file:
                    for category in categories:
                        category_dir = output_file.mkdir(category.name)
                        category_dir.cd()
                        sf_dict[category].Write()
