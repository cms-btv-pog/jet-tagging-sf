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
    is_configured = luigi.BoolParameter(description="Asserts that all parameters and "
        "requirements are set as intended. Used to make sure this task is only executed manually.")

    maximum_error = 0.5 # allowed relative error on scale factors

    def __init__(self, *args, **kwargs):
        super(OptimizeBinning, self).__init__(*args, **kwargs)

    def requires(self):
        return MergeHistograms.req(self, branch=0, version=self.get_version(MergeHistograms),
                _prefer_cli=["version"], iteration=0)

    def store_parts(self):
        return super(OptimizeBinning, self).store_parts() + (self.b_tagger,)

    def output(self):
        return self.wlcg_target("binning.json")

    def run(self):
        if not self.is_configured:
            raise Exception("Cannot run OptimizeBinning task if not run with 'is_configured'")
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

                # determine minimum and maximum bin sizes based on b-tagger working point
                region = category.get_aux("region")
                if self.n_target_bins < 0:
                    # -1 for bin < 0, -1 for number of bin edges vs. bins
                    n_target_bins = len(self.config_inst.get_aux("binning")[region][self.b_tagger]["measurement"]) - 2
                else:
                    n_target_bins = self.n_target_bins

                # figure out best binning to have a reasonable number of events in each bin,
                # as well as a good enough S/B
                signal_hist = hists["signal"]
                bg_hist = hists["bg"]

                # determine rescale factor due to normalization if contamination is shifted
                contamination_factors = self.config_inst.get_aux("contamination_factors")
                contamination_down = contamination_factors["{}_down".format(region)]

                signal_integral = signal_hist.Integral()
                bg_integral = bg_hist.Integral()
                rescale_factor = (signal_integral + bg_integral) / \
                    (contamination_down * bg_integral + signal_integral)

                # start with even split of bins to arrive at *target_bins*
                n_initial_bins = signal_hist.GetNbinsX()
                initial_widths = [signal_hist.GetBinWidth(idx) for idx in range(1, n_initial_bins + 1)]
                bins = np.arange(n_initial_bins)
                merged_bins = np.array_split(bins, n_target_bins)

                starting_widths = [sub_bins[-1] - sub_bins[0] + 1 for sub_bins in merged_bins]

                # helper functions to calculate losses
                def bin_loss(n_signal, err_signal, n_bg, err_bg, region, rescale_factor):
                    # expected error on scale factor, taking into account statistical
                    # and contamination uncertainties

                    # we don't want bins without signal, caught here to prevent ZeroDivision errors
                    if n_signal == 0:
                        return np.inf

                    n_expected = n_signal + n_bg

                    # take into account re-normalization for effect of contamination
                    contamination_error = ((1. - 1./rescale_factor) * n_expected - (1. - contamination_down)*n_bg) / n_signal

                    return np.sqrt(
                        (1. / n_signal)**2 * ( # common factor
                        err_signal**2   # statistical error on signal
                        + n_expected # expected statistical error on data
                        + err_bg**2 # statistical error on contamination
                        )
                        + contamination_error**2 # systematic error on contamination
                    )

                def calculate_losses(widths, signal_hist, bg_hist, rescale_factor):
                    # get bin edges
                    merged_edges = np.cumsum(widths) # one over high edge, for indexing
                    merged_edges = np.insert(merged_edges, 0, 0)

                    edge_positions = np.cumsum(initial_widths)
                    edge_positions = np.insert(edge_positions, 0, 0)
                    bin_edges = edge_positions[merged_edges]

                    # rebin histograms
                    signal_hist_rebinned = signal_hist.Rebin(len(bin_edges) - 1,
                        "signal_rebinned", bin_edges)
                    bg_hist_rebinned = bg_hist.Rebin(len(bin_edges) - 1,
                        "bg_rebinned", bin_edges)

                    # for each bin, calculate expected error on scale factor
                    losses = []
                    base_losses = []
                    for bin_idx in range(1, signal_hist_rebinned.GetNbinsX() + 1):
                        n_signal = signal_hist_rebinned.GetBinContent(bin_idx)
                        n_bg = bg_hist_rebinned.GetBinContent(bin_idx)
                        err_signal = signal_hist_rebinned.GetBinError(bin_idx)
                        err_bg = bg_hist_rebinned.GetBinError(bin_idx)

                        loss = bin_loss(n_signal, err_signal, n_bg, err_bg, region, rescale_factor)
                        base_losses.append(loss) # relative uncertainty
                        losses.append(loss / signal_hist_rebinned.GetBinWidth(bin_idx)) # heuristic to optimize
                    return losses, base_losses

                widths = starting_widths[:]
                best_widths = starting_widths[:]
                best_loss = np.inf
                patience = n_initial_bins
                iter_since_best = 0
                # Optimize binning by calculating a loss value for each bin
                # and increasing the size of the bins with large losses and decreasing
                # the size of those with low losses
                # Break loop if total loss has not improved for *patience* iterations.
                print "Optimizing category {}".format(category.name)
                while True:
                    bin_losses, base_losses = calculate_losses(widths, signal_hist, bg_hist, rescale_factor)
                    loss = sum(bin_losses)

                    if loss < best_loss and max(base_losses) < self.maximum_error:
                        best_widths = widths[:]
                        best_loss = loss
                        iter_since_best = 0
                    else:
                        iter_since_best += 1
                    if iter_since_best > patience:
                        if best_loss == np.inf:
                            raise Exception("Binning optimization failed")
                        break

                    # first make sure that none of the expected errors are larger than a certain value
                    # to reduce likelihood of shifted scale factors going to 0
                    if max(base_losses) > self.maximum_error:
                        sorted_loss_idxs = np.argsort(base_losses)
                    # then optimize heurisitcs
                    else:
                        sorted_loss_idxs = np.argsort(bin_losses)

                    # decrease size of bin with lowest loss that is not already at minimum size
                    for idx in sorted_loss_idxs:
                        if widths[idx] > 1:
                            widths[idx] -= 1
                            break
                    # increase size of bin with highest loss that is not already at maximum size
                    for idx in sorted_loss_idxs[::-1]:
                        if widths[idx] < n_initial_bins:
                            widths[idx] += 1
                            break

                merged_edges = np.cumsum(best_widths)
                bin_edge_values = np.round(merged_edges * signal_hist.GetBinWidth(1), 3)
                results[category.name] = bin_edge_values.tolist()

        print ",\n".join(['"{}": {}'.format(cat_name, binning) for cat_name, binning in results.items()])
        self.output().dump(results, formatter="json", indent=4)
