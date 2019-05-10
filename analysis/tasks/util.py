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
    n_target_bins = luigi.IntParameter(default=-1)
    binning = MergeHistograms.binning # should use even binning here

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

        results = {}

        with inp.load("r") as input_file:
            for category in categories:
                hists = {}

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
                    variable = "jet{}_{}_{}_nominal".format(i_probe_jet, btag_var, region)

                    for process_key in category_dir.GetListOfKeys():
                        process = self.config_inst.get_process(process_key.GetName())
                        if process.is_data or flavor == "inclusive":
                            continue

                        process_dir = category_dir.Get(process.name)
                        hist = process_dir.Get(variable)

                        # separate into signal and contamination
                        if (flavor == "b" and region == "hf") or (flavor == "udsg" and region == "lf"):
                            contrib_type = "signal"
                        else:
                            contrib_type = "bg"

                        if hists.get(contrib_type, None) is None:
                            hists[contrib_type] = hist.Clone()
                        else:
                            hists[contrib_type].Add(hist)

                def bin_loss(n_signal, err_signal, n_bg, err_bg, region):
                    # expected error on scale factor, taking into account statistical
                    # and contamination uncertainties
                    contamination_factors = self.config_inst.get_aux("contamination_factors")
                    contamination_down = contamination_factors["{}_down".format(region)]

                    n_expected = n_signal + n_bg

                    # take into account re-normalization for effect of contamination
                    rescale_factor = n_expected / (contamination_down*n_bg + n_signal)
                    contamination_error = ((1. - 1./rescale_factor) * n_expected - (1. - contamination_down)*n_bg) / n_signal

                    return np.sqrt(
                        (1. / n_signal)**2 * ( # common factor
                        err_signal**2   # statistical error on signal
                        + n_expected # expected statistical error on data
                        + (err_bg**2 + contamination_error**2) # systematic + statistical error on contamination
                    ))

                # determine minimum and maximum bin sizes based on b-tagger working point
                region = category.get_aux("region")
                if self.n_target_bins < 0:
                    # -1 for bin < 0, -1 for number of bin edges vs. bins
                    n_target_bins = len(self.config_inst.get_aux("binning")[region][self.b_tagger]["measurement"]) - 2
                else:
                    n_target_bins = self.n_target_bins

                medium_wp = self.config_inst.get_aux("working_points")[self.b_tagger]["medium"]
                if region == "hf":
                    max_bin_size = medium_wp
                    min_bin_size = (1. - medium_wp) / (2. * n_target_bins)
                elif region == "lf":
                    max_bin_size = 1. - medium_wp
                    min_bin_size = medium_wp / (2. * n_target_bins)
                else:
                    raise ValueError("Unknown region {}.".formt(region))

                # figure out best binning to have a reasonable number of events in each bin,
                # as well as a good enough S/B
                signal_hist = hists["signal"]
                bg_hist = hists["bg"]

                bin_width = signal_hist.GetBinWidth(1)
                signal_values, bg_values = [], []
                signal_errors, bg_errors = [], []

                for bin_idx in range(1, signal_hist.GetNbinsX() + 1):
                    signal_values.append(signal_hist.GetBinContent(bin_idx))
                    signal_errors.append(signal_hist.GetBinError(bin_idx))
                    bg_values.append(bg_hist.GetBinContent(bin_idx))
                    bg_errors.append(bg_hist.GetBinError(bin_idx))

                    if not np.isclose(signal_hist.GetBinWidth(bin_idx), bin_width):
                        raise Exception("Expect even binning to optimize from, but got {}, {}".format(
                            bin_width, signal_hist.GetBinWidth(bin_idx)
                        ))

                # start with even split of bins to arrive at *target_bins*
                n_initial_bins = len(signal_values)
                bins = np.arange(n_initial_bins)
                merged_bins = np.array_split(bins, n_target_bins)
                initial_widths = [sub_bins[-1] - sub_bins[0] + 1 for sub_bins in merged_bins]

                def calculate_losses(widths):
                    # get bin edges
                    merged_high_edges = np.cumsum(widths) # one over high edge, for indexing
                    merged_high_edges[-1] = n_initial_bins # make sure all bins are included
                    merged_low_edges = np.insert(merged_high_edges[:-1], 0, 0)

                    # for each bin, calculate expected error on scale factor
                    losses = []
                    for low_edge, high_edge in zip(merged_low_edges, merged_high_edges):
                        n_signal = sum(signal_values[low_edge:high_edge])
                        err_signal = sum(signal_errors[low_edge:high_edge])
                        n_bg = sum(bg_values[low_edge:high_edge])
                        err_bg = sum(bg_errors[low_edge:high_edge])

                        loss = bin_loss(n_signal, err_signal, n_bg, err_bg, region)
                        # penalty term to enforce more even binning
                        losses.append(loss / float(high_edge - low_edge))

                    return losses

                # translate bin size limits to limits on the number of bins
                min_merged_bins = min_bin_size * n_initial_bins
                max_merged_bins = max_bin_size * n_initial_bins

                widths = initial_widths[:]
                best_widths = initial_widths[:]
                best_loss = np.inf
                patience = 10
                iter_since_best = 0
                # Optimize binning by calculating a loss value for each bin
                # and increasing the size of the bins with large losses and decreasing
                # the size of those with low losses
                # Break loop if total loss has not improved for *patience* iterations.
                while True:
                    bin_losses = calculate_losses(widths)
                    loss = sum(bin_losses)

                    if loss < best_loss:
                        best_widths = widths[:]
                        best_loss = loss
                        iter_since_best = 0
                    else:
                        iter_since_best += 1
                    if iter_since_best > patience:
                        break

                    sorted_loss_idxs = np.argsort(bin_losses)

                    # decrease size of bin with lowest loss that is not already at minimum size
                    for idx in sorted_loss_idxs:
                        if widths[idx] > min_merged_bins:
                            widths[idx] -= 1
                            break
                    # increase size of bin with highest loss that is not already at maximum size
                    for idx in sorted_loss_idxs[::-1]:
                        if widths[idx] < max_merged_bins:
                            widths[idx] += 1
                            break

                merged_edes = np.cumsum(best_widths)
                bin_edge_values = merged_edes * bin_width
                results[category.name] = [round(edge, 3) for edge in bin_edge_values]
                print [round(edge, 3) for edge in bin_edge_values]

        self.output().dump(results, formatter="json", indent=4)
