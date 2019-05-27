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

from analysis.tasks.base import AnalysisTask
from analysis.root import ROOTPlot
from law.target.local import LocalDirectoryTarget

dirname = os.path.abspath(os.path.dirname(__file__))

class PlotFromCSV(AnalysisTask):
    shift = luigi.Parameter()
    flavor = luigi.ChoiceParameter(choices=["lf", "c", "hf"])
    csv_file = os.path.join(dirname, "DeepCSV_94XSF_V3_B_F.csv")
    compare_file = os.path.join(dirname, "DeepCSV_94XSF_V4_B_F.csv")

    def output(self):
        return self.local_target("plots.tgz")

    def run(self):
        ROOT.PyConfig.IgnoreCommandLineOptions = True
        ROOT.gROOT.SetBatch()

        ROOT.gSystem.Load('libCondFormatsBTauObjects')
        ROOT.gSystem.Load('libCondToolsBTau')

        local_tmp = LocalDirectoryTarget(is_tmp=True)
        local_tmp.touch()

        v_sys = getattr(ROOT, 'vector<string>')()
        v_sys.push_back("up_" + self.shift)
        v_sys.push_back("down_" + self.shift)

        readers = {}
        for input_file, id in zip([self.csv_file, self.compare_file], ["first", "compare"]):
            calib = ROOT.BTagCalibration("csvv1", input_file)
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
            readers[id] = reader

        flavor_ids = self.config_inst.get_aux("flavor_ids")
        binning = self.config_inst.get_aux("binning")[self.flavor]

        pt_binning = [(start, end) for start, end in zip(binning["pt"][:-1], binning["pt"][1:])]
        eta_binning = [(start, end) for start, end in zip(binning["abs(eta)"][:-1], binning["abs(eta)"][1:])]

        for pt_idx, pt_range in enumerate(pt_binning):
            for eta_idx, eta_range in enumerate(eta_binning):
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_title("pt: %s to %s, eta: %.1f to %.1f" % (pt_range + eta_range))

                pt_val = np.mean(pt_range)
                eta_val = np.mean(eta_range)

                for sys_type in ["central", "up_" + self.shift, "down_" + self.shift]:
                    for id, reader in readers.items():
                        x_values = np.linspace(0., 1., 10000)
                        y_values = []
                        for csv_value in x_values:
                            sf = reader.eval_auto_bounds(
                                sys_type,      # systematic (here also 'up'/'down' possible)
                                flavor_ids[self.flavor],              # jet flavor
                                eta_val,            # absolute value of eta
                                pt_val,             # pt
                                csv_value
                            )
                            y_values.append(sf)
                        ax.plot(x_values, y_values, label="{}, {}".format(id, sys_type))

                ax.legend()
                fig.savefig(os.path.join(local_tmp.path, "SF_%s_%s_Pt%sto%s_eta%.1fTo%.1f.pdf" % ((self.flavor, self.shift) + pt_range + eta_range)))

        with self.output().localize("w") as tmp:
            with tarfile.open(tmp.path, "w:gz") as tar:
                for plot_file in os.listdir(local_tmp.path):
                    tar.add(os.path.join(local_tmp.path, plot_file), arcname=plot_file)


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
