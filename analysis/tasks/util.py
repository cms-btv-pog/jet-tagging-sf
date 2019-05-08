# -*- coding: utf-8 -*-

import law
import luigi
import numpy as np

from collections import defaultdict
from scipy.optimize import minimize

from analysis.config.jet_tagging_sf import get_category
from analysis.tasks.base import AnalysisTask
from analysis.tasks.hists import MergeHistograms

class OptimizeBinning(AnalysisTask):

    category_tags = MergeHistograms.category_tags

    b_tagger = MergeHistograms.b_tagger
    n_target_bins = luigi.IntParameter(default=12)
    n_bins = MergeHistograms.n_bins

    min_bin_size = 0.02
    max_bin_size = 0.5

    def requires(self):
        return MergeHistograms.req(self, branch=0, version=self.get_version(MergeHistograms),
                _prefer_cli=["version"], iteration=0)

    def store_parts(self):
        return super(OptimizeBinning, self).store_parts() + (self.b_tagger,)

    def output(self):
        return self.wlcg_target("binning.json")

    def run(self):
        import ROOT

        inp = self.input()

        categories = []
        for category, _, _ in self.config_inst.walk_categories():
            if category.has_tag(("merged", self.b_tagger), mode=all) and category.get_aux("phase_space") == "measure":
                if len(self.category_tags) > 0 and not category.has_tag(self.category_tags, mode=any):
                    continue
                categories.append(category)

        with inp.load("r") as input_file:
            # region -> signal/contamination -> histogram
            hists = defaultdict(dict)

            for category in categories:
                for leaf_cat, _, children in category.walk_categories():
                    # we are only interested in leaves
                    if children:
                        continue

                    category_dir = input_file.Get(leaf_cat.name)

                    flavor = leaf_cat.get_aux("flavor")
                    region = leaf_cat.get_aux("region")
                    i_probe_jet = leaf_cat.get_aux("i_probe_jet")

                    # Add up b-tag discriminants for probe jets in each region, split by flavor
                    btag_var = self.config_inst.get_aux("btaggers")[self.b_tagger]["variable"]
                    variable = "jet{}_{}_{}_nominal".format(i_probe_jet, btag_var, region.upper())

                    for process_key in category_dir.GetListOfKeys():
                        process = self.config_inst.get_process(process_key.GetName())
                        if process.is_data:
                            continue

                        process_dir = category_dir.Get(process.name)
                        hist = process_dir.Get(variable)

                        # separate into signal and contamination
                        if (flavor == "b" and region == "hf") or (flavor == "udsg" and region == "lf"):
                            contrib_type = "signal"
                        else:
                            contrib_type = "bg"

                        if hists[region].get(contrib_type, None) is None:
                            hists[region][contrib_type] = hist.Clone()
                        else:
                            hists[region][contrib_type].Add(hist)

            def bin_loss(n_signal, err_signal, n_bg, err_bg, region):
                # expected error on scale factor, taking into account statistical
                # and contamination uncertainties
                contamination_factors = self.config_inst.get_aux("contamination_factors")
                contamination_shift = 1. - contamination_factors["{}_down".format(region)]
                return np.sqrt(
                    (1. / n_signal)**2 * ( # common factor
                    err_signal**2   # statistical error on signal
                    + (n_signal + n_bg) # expected statistical error on data
                    + (err_bg**2 + (contamination_shift*n_bg)**2) # systematic + statistical error on contamination
                ))

            # figure out best binning to have a reasonable number of events in each bin,
            # as well as a good enough S/B
            for region in ["hf", "lf"]:
                signal_hist = hists[region]["signal"]
                bg_hist = hists[region]["bg"]

                bin_edges = []
                signal_values, bg_values = [], []
                signal_errors, bg_errors = [], []

                for bin_idx in range(1, signal_hist.GetNbinsX() + 1): # TODO: Only bins larger than zero
                    signal_values.append(signal_hist.GetBinContent(bin_idx))
                    signal_errors.append(signal_hist.GetBinError(bin_idx))
                    bg_values.append(bg_hist.GetBinContent(bin_idx))
                    bg_errors.append(bg_hist.GetBinError(bin_idx))

                    bin_edges.append(bg_hist.GetBinLowEdge(bin_idx))

                # start with even split of bins to arrive at *target_bins*
                bins = np.arange(self.n_bins)
                merged_bins = np.array_split(bins, self.n_target_bins)
                initial_widths = [sub_bins[-1] - sub_bins[0] + 1 for sub_bins in merged_bins]

                def merge_loss(widths):
                    # get bin edges
                    merged_high_edges = np.cumsum(widths) # one over high edge, for indexing
                    merged_low_edges = np.insert(merged_high_edges[:-1], 0, 0)

                    # for each bin, calculate expected error on scale factor
                    losses = []
                    for low_edge, high_edge in zip(merged_low_edges, merged_high_edges):
                        n_signal = sum(signal_values[low_edge:high_edge])
                        err_signal = sum(signal_errors[low_edge:high_edge])
                        n_bg = sum(bg_values[low_edge:high_edge])
                        err_bg = sum(bg_errors[low_edge:high_edge])

                        loss = bin_loss(n_signal, err_signal, n_bg, err_bg, region)
                        losses.append(loss)

                    # metric: worst bin score
                    return max(losses)

                # perform minimization, taking into acount minimum and maximum bin width
                # constrain sum of widths to 1
                constraints = {
                    "type": "eq",
                    "fun": lambda x: sum(x) - 1,
                }
                result = minimize(merge_loss, initial_widths, constraints=constraints,
                    bounds=((self.min_bin_size, self.max_bin_size),) * self.n_target_bins
                )



        self.output().dump(binning, formatter="json", indent=4)
