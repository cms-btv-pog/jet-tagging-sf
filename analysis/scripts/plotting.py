# -*- coding: utf-8 -*-

import ROOT

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"figure.autolayout": True})
plt.style.use("ggplot")

import numpy as np
import itertools
import os
import tarfile

import luigi
import law

from collections import OrderedDict
from analysis.tasks.base import AnalysisTask
from analysis.root import ROOTPlot
from law.target.local import LocalDirectoryTarget

dirname = os.path.abspath(os.path.dirname(__file__))

class PlotFromCSV(AnalysisTask):
    shift = luigi.Parameter(default="ALL")
    flavor = luigi.ChoiceParameter(choices=["lf", "c", "hf"])
    csv_file = os.path.join(dirname, "sf_2018_prod11.csv")  # new
    compare_file = os.path.join(dirname, "sf_2018_prod7.csv")  # old
    norm_to_nominal = luigi.BoolParameter()

    root_hf_file = os.path.join(dirname, "scale_factors_deepjet_hf_binned.root")
    root_lf_file = os.path.join(dirname, "scale_factors_deepjet_lf_binned.root")

    def output(self):
        return self.local_target("plots_{}.tgz".format(self.shift))

    def run(self):
        ROOT.PyConfig.IgnoreCommandLineOptions = True
        ROOT.gROOT.SetBatch()

        ROOT.gSystem.Load('libCondFormatsBTauObjects')
        ROOT.gSystem.Load('libCondToolsBTau')

        local_tmp = LocalDirectoryTarget(is_tmp=True)
        local_tmp.touch()

        jes_sources = self.config_inst.get_aux("jes_sources")
        shifts = []
        if self.shift == "ALL":
            if self.flavor == "c":
                shifts.extend(["cferr1", "cferr2"])
            else:
                shifts.extend(["jes{}".format(jes_source) for jes_source in jes_sources if jes_source != "Total"])
                shifts.extend(["{}{}".format(region, type) for region, type in
                    itertools.product(["lf", "hf"], ["", "stats1", "stats2"])])
        elif self.shift == "NONE":
            shifts = []
        else:
            shifts = [self.shift]

        v_sys = getattr(ROOT, 'vector<string>')()
        for shift in shifts:
            v_sys.push_back("up_" + shift)
            v_sys.push_back("down_" + shift)

        flavor_ids = self.config_inst.get_aux("flavor_ids")
        binning = self.config_inst.get_aux("binning")[self.flavor]

        pt_binning = [(start, end) for start, end in zip(binning["pt"][:-1], binning["pt"][1:])]
        eta_binning = [(start, end) for start, end in zip(binning["abs(eta)"][:-1], binning["abs(eta)"][1:])]

        figures = {}
        if self.compare_file is None:
            csv_files, descriptions = [self.csv_file], ["csv"]
        else:
            csv_files, descriptions = [self.csv_file, self.compare_file], ["new", "old"]
        for input_file, id in zip(csv_files, descriptions):
            # create calibration reader
            calib = ROOT.BTagCalibration("csv_{}".format(id), input_file)
            reader = ROOT.BTagCalibrationReader(
                3,              # 0 is for loose op, 1: medium, 2: tight, 3: discr. reshaping
                "central",      # central systematic type
                v_sys,          # vector of other sys. types
            )
            for jetFlavor in [0, 1, 2]:
                reader.load(
                    calib,
                    jetFlavor,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG
                    "iterativefit"      # measurement type
                )

            for pt_idx, pt_range in enumerate(pt_binning):
                for eta_idx, eta_range in enumerate(eta_binning):
                    key = pt_range + eta_range
                    if key not in figures:
                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        ax.set_title("pt: %s to %s, eta: %.1f to %.1f" % (pt_range + eta_range))
                    else:
                        fig, ax = figures[key]

                    if pt_range[1] == np.inf:
                        pt_val = pt_range[0] + 1
                    else:
                        pt_val = np.mean(pt_range)
                    eta_val = np.mean(eta_range)

                    def get_values(csv_reader, sys_type):
                        x_values = np.linspace(-0.1, 1., 10000)
                        y_values = []
                        for csv_value in x_values:
                            sf = csv_reader.eval_auto_bounds(
                                sys_type,      # systematic (here also 'up'/'down' possible)
                                flavor_ids[self.flavor],              # jet flavor
                                eta_val,            # absolute value of eta
                                pt_val,             # pt
                                csv_value
                            )
                            y_values.append(sf)
                        return np.array(x_values), np.array(y_values)

                    x_values, nominal_values = get_values(reader, "central")
                    if not self.norm_to_nominal:
                        ax.plot(x_values, nominal_values, label="{}, {}".format(id, "nominal"))

                    if self.shift != "NONE":
                        total_errors_up = np.zeros(nominal_values.shape)
                        total_errors_down = np.zeros(nominal_values.shape)
                        for shift in shifts:
                            _, up_values = get_values(reader, "up_" + shift)
                            _, down_values = get_values(reader, "down_" + shift)

                            if len(shifts) > 1: # build envelope
                                diff_up = up_values - nominal_values
                                diff_down = down_values - nominal_values

                                # shift with effect in up/down direction
                                errors_up = np.max([diff_up, diff_down, np.zeros(nominal_values.shape)], axis=0)
                                errors_down = np.min([diff_up, diff_down, np.zeros(nominal_values.shape)], axis=0)

                                # add in quadrature
                                total_errors_up += errors_up**2
                                total_errors_down += errors_down**2
                        total_errors_up = total_errors_up**0.5
                        total_errors_down = total_errors_down**0.5

                        if len(shifts) > 1:
                            up_values = nominal_values + total_errors_up
                            down_values = nominal_values - total_errors_down
                        if self.norm_to_nominal:
                            up_values /= nominal_values
                            down_values /= nominal_values
                        ax.plot(x_values, up_values, label="{}, {}".format(id, "up_" + self.shift))
                        ax.plot(x_values, down_values, label="{}, {}".format(id, "down" + self.shift))

                    if self.compare_file is None:
                        if self.flavor in ["c", "hf"]:
                            input_file = self.root_hf_file
                        elif self.flavor == "lf":
                            input_file = self.root_lf_file
                        else:
                            raise Exception("No .root file for c flavor SFs.")

                        root_file = ROOT.TFile.Open(input_file)

                        func_name = "csv_ratio_Pt{}_Eta{}_final".format(pt_idx, eta_idx)
                        if self.flavor == "c":
                            func_name = "c_" + func_name

                        func = root_file.Get(func_name)

                        y_values = []
                        x_values = np.linspace(-0.1, 1., 10000)
                        for csv_value in x_values:
                            y_val = func.Eval(csv_value)
                            y_values.append(y_val)

                        ax.plot(x_values, y_values, label="{}, {}".format(".root", "nominal"))

                    figures[key] = (fig, ax)

            del reader
            del calib

        for key, (fig, ax) in figures.items():
            ax.legend(loc="lower right")
            ax.set_ylim(0., 2.)
            fig.savefig(os.path.join(local_tmp.path, "SF_%s_%s_Pt%sTo%s_eta%.1fTo%.1f.pdf" % ((self.flavor, self.shift) + key)))

        with self.output().localize("w") as tmp:
            with tarfile.open(tmp.path, "w:gz") as tar:
                for plot_file in os.listdir(local_tmp.path):
                    tar.add(os.path.join(local_tmp.path, plot_file), arcname=plot_file)


class PlotShiftsFromCSV(AnalysisTask, law.WrapperTask):
    shifts = law.CSVParameter(default=[], description="shifts to require")
    skip_shifts = law.CSVParameter(default=[], description="shifts to skip, supports patterns")

    flavor = PlotFromCSV.flavor
    norm_to_nominal = PlotFromCSV.norm_to_nominal

    wrapped_task = PlotFromCSV

    def __init__(self, *args, **kwargs):
        super(PlotShiftsFromCSV, self).__init__(*args, **kwargs)

        jes_sources = self.config_inst.get_aux("jes_sources")

        if not self.shifts:
            self.shifts = []
            if self.flavor == "c":
                self.shifts = ["cferr1", "cferr2"]
            else:
                self.shifts.extend(["jes{}".format(jes_source) for jes_source in jes_sources if jes_source != "Total"])
                self.shifts.extend(["{}{}".format(region, type) for region, type in
                    itertools.product(["lf", "hf"], ["", "stats1", "stats2"])])
        if self.skip_shifts:
            filter_fn = lambda d: not law.util.multi_match(d, self.skip_shifts)
            self.shifts = filter(filter_fn, self.shifts)

    def requires(self):
        def req(shift):
            return self.wrapped_task.req(self, shift=shift)

        return OrderedDict([(shift, req(shift)) for shift in self.shifts])


class PlotFromRoot(AnalysisTask):

    hf_file = "/user/rath/Deepcsv_rwt_fit_hf_v2_final_2018_2_12test.root"
    lf_file = "/user/rath/Deepcsv_rwt_fit_lf_v2_final_2018_2_12test.root"

    flavor = luigi.ChoiceParameter(choices=["hf", "lf"])
    norm_to_nominal = luigi.BoolParameter()
    shift = luigi.Parameter(default="ALL")

    def output(self):
        return self.local_target("plots.tgz")

    def run(self):
        ROOT.PyConfig.IgnoreCommandLineOptions = True
        ROOT.gROOT.SetBatch()

        local_tmp = LocalDirectoryTarget(is_tmp=True)
        local_tmp.touch()

        if self.flavor == "hf":
            input_file = self.hf_file
        else:
            input_file = self.lf_file

        root_file = ROOT.TFile.Open(input_file)

        hist_tpl = "h_csv_ratio_Pt{}_Eta{}_final"

        binning = self.config_inst.get_aux("binning")[self.flavor]
        n_pt_categories = len(binning["pt"]) - 1
        n_eta_categories = len(binning["abs(eta)"]) - 1
        hist_names = [hist_tpl.format(pt_idx, eta_idx)
            for pt_idx in range(n_pt_categories) for eta_idx in range(n_eta_categories)]

        ##
        ## Plot histograms with total systematic envelope
        ##

        for hist_name in hist_names:
            nominal_hist = root_file.Get(hist_name)

            errors_up = []
            errors_down = []
            # collect all shifts from file
            if self.shift == "ALL":
                shifts = [key.GetName().split("_")[-1] for key in root_file.GetListOfKeys()
                    if key.GetName().startswith(hist_name)]
                shifts = set([shift[:-2] for shift in shifts if shift[-2:] == "Up"])
            else:
                shifts = [self.shift]
            for shift_idx, shift in enumerate(shifts):
                hist_name_up = hist_name + "_" + shift + "Up"
                hist_name_down = hist_name + "_" + shift + "Down"
                hist_up = root_file.Get(hist_name_up)
                hist_down = root_file.Get(hist_name_down)

                for bin_idx in range(1, nominal_hist.GetNbinsX() + 1):
                    nominal_value = nominal_hist.GetBinContent(bin_idx)

                    # combine all shifts that have an effect in the same direction
                    # effect from <shift>_up/done systematics
                    diff_up = hist_up.GetBinContent(bin_idx) - nominal_value
                    diff_down = hist_down.GetBinContent(bin_idx) - nominal_value


                    # detect systematics where up/down shift direction is the same
                    #if diff_up * diff_down > 0:
                    #    print "One sided shift: {}, {}".format(shift, category)

                    # if multiple shifts, build envelope
                    if len(shifts) != 1:
                        # shift with effect in up/down direction
                        error_up = max([diff_up, diff_down, 0])
                        error_down = min([diff_up, diff_down, 0])

                        # add in quadrature
                        if shift_idx == 0:
                            errors_up.append(error_up**2)
                            errors_down.append(error_down**2)
                        else:
                            errors_up[bin_idx - 1] += error_up**2
                            errors_down[bin_idx - 1] += error_down**2
                    else:
                        errors_up.append(diff_up)
                        errors_down.append(-diff_down) # is subtracted later
            # multiple shifts have been added quadratically, take square root
            if len(shifts) != 1:
                errors_up = np.sqrt(errors_up)
                errors_down = np.sqrt(errors_down)

            # build shifted histograms
            combined_hist_up = nominal_hist.Clone()
            combined_hist_down = nominal_hist.Clone()

            for bin_idx in range(1, nominal_hist.GetNbinsX() + 1):
                combined_hist_up.SetBinContent(bin_idx, combined_hist_up.GetBinContent(bin_idx)
                    + errors_up[bin_idx - 1])
                combined_hist_down.SetBinContent(bin_idx, combined_hist_down.GetBinContent(bin_idx)
                    - errors_down[bin_idx - 1])

            if self.norm_to_nominal:
                combined_hist_up.Divide(nominal_hist)
                combined_hist_down.Divide(nominal_hist)

            plot = ROOTPlot(hist_name, hist_name)
            plot.create_pads()
            plot.cd(0, 0)
            plot.draw({"nominal": nominal_hist}, line_color=1)
            plot.draw({"up": combined_hist_up}, line_color=2)
            plot.draw({"down": combined_hist_down}, line_color=4)

            plot.save(os.path.join(local_tmp.path, "{}.pdf".format(hist_name)))

        ##
        ## Check scale factors, uncertainties, and fits
        ##

        for i, (pt_idx, eta_idx) in enumerate(itertools.product(range(n_pt_categories),
            range(n_eta_categories))):
            data_hist = root_file.Get("h_csv_Data_Pt{}_Eta{}".format(pt_idx, eta_idx))
            if self.flavor == "hf":
                signal_base = "h_csv_MC_bjets"
                bg_base = "h_csv_MC_nonbjets"
            else:
                signal_base = "h_csv_MC_nonbjets" # actually lf
                bg_base = "h_csv_MC_bjets" # actually b + c

            signal_hist = root_file.Get("{}_Pt{}_Eta{}".format(signal_base, pt_idx, eta_idx))
            bg_hist = root_file.Get("{}_Pt{}_Eta{}".format(bg_base, pt_idx, eta_idx))

            #

        with self.output().localize("w") as tmp:
            with tarfile.open(tmp.path, "w:gz") as tar:
                for plot_file in os.listdir(local_tmp.path):
                    tar.add(os.path.join(local_tmp.path, plot_file), arcname=plot_file)
