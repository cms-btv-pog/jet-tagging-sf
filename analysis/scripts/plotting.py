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

from analysis.tasks.base import AnalysisTask
from law.target.local import LocalDirectoryTarget

dirname = os.path.abspath(os.path.dirname(__file__))

class PlotFromCSV(AnalysisTask):

    inputFile = os.path.join(dirname, "DeepCSV_94XSF_V3_B_F.csv")

    def output(self):
        return self.local_target("plots.tgz")

    def run(self):
        ROOT.PyConfig.IgnoreCommandLineOptions = True
        ROOT.gROOT.SetBatch()

        ROOT.gSystem.Load('libCondFormatsBTauObjects')
        ROOT.gSystem.Load('libCondToolsBTau')

        outp = self.output()
        local_tmp = LocalDirectoryTarget(is_tmp=True)
        local_tmp.touch()

        calib = ROOT.BTagCalibration("csvv1", self.inputFile)

        v_sys = getattr(ROOT, 'vector<string>')()
        v_sys.push_back('up_lf')
        v_sys.push_back('down_lf')

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

        flavors = {"LF": 2, "HF": 0}
        binning = self.config_inst.get_aux("binning")
        for flavorName, flavorID in flavors.items():
            flavor_bins = binning[flavorName]
            csv_bins = flavor_bins["deepcsv"]

            pt_binning = [(start, end) for start, end in zip(flavor_bins["pt"][:-1], flavor_bins["pt"][1:])]
            eta_binning = [(start, end) for start, end in zip(flavor_bins["abs(eta)"][:-1], flavor_bins["abs(eta)"][1:])]
       	    csv_binning = [(start, end) for start, end in zip(csv_bins[:-1], csv_bins[1:])]

            csv_min = min(csv_bins)
       	    csv_max = max(csv_bins)

            for eta_range, pt_range in itertools.product(eta_binning, pt_binning):
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_title("pt: %s to %s, eta: %.1f to %.1f" % (pt_range + eta_range))
                
                eta_val = sum(eta_range)/2.
                pt_val = sum(pt_range)/2.

                x_values = np.linspace(csv_min, csv_max, 10000)
                y_values = []
                for csv_value in x_values:
                    sf = reader.eval_auto_bounds(
                        'central',      # systematic (here also 'up'/'down' possible)
                        flavorID,              # jet flavor
                        eta_val,            # absolute value of eta
                        pt_val,             # pt
                        csv_value
                    )
                    y_values.append(sf)

                ax.plot(x_values, y_values)
                fig.savefig(os.path.join(local_tmp.path, "SF_%s_Pt%sto%s_eta%.1fTo%.1f.pdf" % ((flavorName,) + pt_range + eta_range)))
                
        with outp.localize("w") as tmp:
            with tarfile.open(tmp.path, "w:gz") as tar:
                for plot_file in os.listdir(local_tmp.path):
                    tar.add(os.path.join(local_tmp.path, plot_file), arcname=plot_file)
