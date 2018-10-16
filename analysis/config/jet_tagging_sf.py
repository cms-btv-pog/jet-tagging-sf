# -*- coding: utf-8 -*-
# flake8: noqa

"""
Definition of the analysis for extracting jet tagging scale factors
as well as its configurations for different campaigns.
"""


import re

import numpy as np
import order as od
import scinum as sn
import six
from order.util import join_root_selection

from analysis.config.constants import Z_MASS
from analysis.config.processes import process_data_ee, process_data_emu, process_data_mumu, \
    process_tt_dl, process_dy_lep, process_st_tW, process_WW_sl


# define the analysis
analysis = od.Analysis("jet_tagging_sf", 1)


# setup the config for ICHEP 2018
from analysis.config.campaign_ICHEP18 import campaign as campaign_ICHEP18
config_ICHEP18 = cfg = analysis.add_config(campaign=campaign_ICHEP18)


# link processes
cfg.add_process(process_data_ee)
cfg.add_process(process_data_emu)
cfg.add_process(process_data_mumu)
cfg.add_process(process_tt_dl)
cfg.add_process(process_dy_lep)
cfg.add_process(process_st_tW)
cfg.add_process(process_WW_sl)

# add datasets
dataset_names = [
    "data_B_ee", "data_C_ee", "data_D_ee", "data_E_ee", "data_F_ee",
    "data_B_emu", "data_C_emu", "data_D_emu", "data_E_emu", "data_F_emu",
    "data_B_mumu", "data_C_mumu", "data_D_mumu", "data_E_mumu", "data_F_mumu",
    "tt_dl",
    #"dy_lep_4To50_Ht70To100",
    #"dy_lep_4To50_Ht100To200", "dy_lep_4To50_Ht200To400",
    #"dy_lep_4To50_Ht400To600", "dy_lep_4To50_Ht600ToInf",
    "dy_lep_10To50",
    #"dy_lep_50ToInf",
    "dy_lep_50ToInf_Ht70To100", "dy_lep_50ToInf_Ht100To200",
    "dy_lep_50ToInf_Ht200To400", "dy_lep_50ToInf_Ht400To600",
    "dy_lep_50ToInf_Ht600To800", "dy_lep_50ToInf_Ht800To1200", "dy_lep_50ToInf_Ht1200To2500",
    "dy_lep_50ToInf_Ht2500ToInf",
    "st_tW_t", "st_tW_tbar",
    "WW_sl",
]
for dataset_name in dataset_names:
    dataset = campaign_ICHEP18.get_dataset(dataset_name)
    cfg.add_dataset(dataset)

# define channels
ch_ee = cfg.add_channel("ee", 1)
ch_emu = cfg.add_channel("emu", 2)
ch_mumu = cfg.add_channel("mumu", 3)

# store channels per real dataset
cfg.set_aux("dataset_channels", {
    dataset: cfg.get_channel(dataset.name.split("_")[-1])
    for dataset in cfg.datasets.values()
    if dataset.is_data
})

# store DeepCSV working points
cfg.set_aux("working_points", {
    "deepcsv": {
        "loose": 0.1522,
        "medium": 0.4941,
        "tight": 0.8001,
    }
})
cfg.set_aux("btagger", {
    "name": "deepcsv",
    "variable": "deepcsv_bcomb",
})

# store binning information
cfg.set_aux("binning", {
    "LF": {
        "pt": [20, 30, 40, 60, np.inf],
        "abs(eta)": [0., 0.8, 1.6, 2.4],
        "deepcsv": {
            "plotting": [
                -2.01, 0.0, 0.0254, 0.0508, 0.0762, 0.1016, 0.127, 0.1522, 0.2205, 0.2889, 0.3573,
                0.4257, 0.4941, 0.5961, 0.6981, 0.8001, 0.835, 0.87, 0.905, 0.94, 0.975, 1.01,
            ],
            "measurement": [
                -2.01, 0.0, 0.0254, 0.0508, 0.0762, 0.1016, 0.127, 0.1522, 0.2205, 0.2889, 0.3573,
                0.4257, 0.4941, 1.01,
            ],
        },
    },
    "HF": {
        "pt": [20, 30, 50, 70, 100, np.inf],
        "abs(eta)": [0., 2.4],
        "deepcsv": {
            "plotting": [
                -2.01, 0.0, 0.0254, 0.0508, 0.0762, 0.1016, 0.127, 0.1522, 0.2205, 0.2889, 0.3573,
                0.4257, 0.4941, 0.5553, 0.6165, 0.6777, 0.7389, 0.8001, 0.842, 0.884, 0.926, 0.968,
                1.01,
            ],
            "measurement": [
                -2.01, 0.0, 0.1522, 0.2205, 0.2889, 0.3573, 0.4257, 0.4941, 0.5553, 0.6165, 0.6777,
                0.7389, 0.8001, 0.842, 0.884, 0.926, 0.968, 1.01
            ],
        },
    }
})

# define nested categories (analysis phase space -> hf/lf region -> flavor -> pt bin -> eta bin)
def get_phasespace_info():
    return [
        ("measure", join_root_selection(["n_jets == 2", "mll > 12", "dr_ll > 0.2"])),
        ("closure", join_root_selection(["n_jets >= 2", "mll > 12", "dr_ll > 0.2"])),
    ]

def get_region_info(idx, channel, et_miss=30., z_window=10., add_btag_cut=True):
    btagger = cfg.get_aux("btagger")
    btag_name = btagger["name"]
    btag_variable = cfg.get_variable("jet{}_{}".format(idx, btagger["variable"]))

    csv_tight = cfg.get_aux("working_points")[btag_name]["tight"]
    csv_loose = cfg.get_aux("working_points")[btag_name]["loose"]

    hf_cuts, lf_cuts = [], []
    if add_btag_cut:
        # jet tagging requirement
        hf_cuts.append("({}) > {}".format(btag_variable.expression, csv_tight))
        lf_cuts.append("({}) < {}".format(btag_variable.expression, csv_loose))

    # the following cuts do not apply for emu
    if channel != "emu":
        # ET-miss requirement
        et_miss_expr = "(met_px**2 + met_py**2)**0.5"
        hf_cuts.append("{} > {}".format(et_miss_expr, et_miss))
        lf_cuts.append("{} < {}".format(et_miss_expr, et_miss))

        # z-mass window
        hf_cuts.append("abs(mll - {}) > {}".format(Z_MASS.nominal, z_window))
        lf_cuts.append("abs(mll - {}) < {}".format(Z_MASS.nominal, z_window))

        # z peak diamond
        lf_cuts.append("pass_z_mask == 0")

    return [
        ("HF", join_root_selection(hf_cuts)),
        ("LF", join_root_selection(lf_cuts)),
    ]

def get_flavor_info(idx):
    return [
        ("b", "abs(jet%d_flavor) == 5" % idx),
        ("c", "abs(jet%d_flavor) == 4" % idx),
        ("udsg", "abs(jet%d_flavor) != 5 && abs(jet%d_flavor) != 4" % (idx, idx)),
        ("inclusive", "1 * 1"),
    ]

def binning_to_selection(binning, variable):
    def to_string(value):
        if np.isinf(value):
            return "Inf"
        elif value == 0:
            return "0"
        else:
            return str(value).replace(".", "p")

    selections = []
    for left_edge, right_edge in zip(binning[:-1], binning[1:]):
        name = "{}To{}".format(to_string(left_edge), to_string(right_edge))
        cuts = []
        if left_edge != 0:
            cuts.append("{} > {}".format(variable, left_edge))
        if not np.isinf(right_edge):
            cuts.append("{} <= {}".format(variable, right_edge))
        edges = [left_edge, right_edge]
        selections.append((name, join_root_selection(cuts), edges))
    return selections

def get_axis_info(idx, axis_var, fmt=None):
    if fmt is None:
        fmt = "jet{}_" + axis_var

    binning = cfg.get_aux("binning")
    hf_bins = binning["HF"][axis_var]
    lf_bins = binning["LF"][axis_var]
    variable = fmt.format(idx)
    return {
        "HF": binning_to_selection(hf_bins, variable),
        "LF": binning_to_selection(lf_bins, variable),
    }

def get_category(pt, eta, region, phase_space="measure"):
    matches = []
    for category in cfg.categories:
        cat_phasespace = category.get_aux("phase_space", None)
        if not cat_phasespace == phase_space:
            continue

        cat_region = category.get_aux("region", None)
        if not region == cat_region:
            continue

        cat_pt_range = category.get_aux("pt", (0., 0.))
        if not (cat_pt_range[0] < pt <= cat_pt_range[1]):
            continue

        cat_eta_range = category.get_aux("eta", (0., 0.))
        if not (cat_eta_range[0] < eta <= cat_eta_range[1]):
            continue
        matches.append(category)
    # only return the category if the matching is unambiguos
    if len(matches) == 1:
        return matches[0]
    else:
        raise ValueError("Expected one single matching category, but got {}".format(matches))

# variables
#cfg.add_variable(
#    name="dr_ll",
#    expression="dr_ll",
#    binning=(25, 0., 5.,),
#    x_title="dR(ll)",
#)
cfg.add_variable(
    name="n_jets",
    expression="n_jets",
    binning=(10, 0., 10.,),
    x_title="N(jets)",
)
for lep_idx in xrange(1, 3):
    cfg.add_variable(
        name="lep{}_pt".format(lep_idx),
        expression="(lep{}_px**2 + lep{}_py**2)**0.5".format(lep_idx, lep_idx),
        binning=(25, 0., 500.,),
        unit="GeV",
        x_title="Lep_{} p_{{T}}".format(lep_idx),
    )

for jet_idx in xrange(1, 5):
    cfg.add_variable(
        name="jet{}_pt".format(jet_idx),
        expression="jet{}_pt".format(jet_idx),
        binning=(25, 0., 500.,),
        unit="GeV",
        x_title="Jet_{} p_{{T}}".format(jet_idx),
    )
    for region in [None, "HF", "LF"]:
        if not region:
            binning = (25, 0., 1.)
            tags = set()
            postfix = ""
        elif region == "HF":
            binning = cfg.get_aux("binning")["HF"]["deepcsv"]["plotting"]
            tags = {"skip_LF"}
            postfix = "_HF"
        elif region == "LF":
            binning = cfg.get_aux("binning")["LF"]["deepcsv"]["plotting"]
            tags = {"skip_HF"}
            postfix = "_LF"

        cfg.add_variable(
            name="jet{}_deepcsv_b{}".format(jet_idx, postfix),
            expression="jet{}_deepcsv_b".format(jet_idx),
            binning=binning,
            x_title="Jet_{} prob_{{b}}".format(jet_idx),
            tags=tags,
        )
        cfg.add_variable(
            name="jet{}_deepcsv_bb{}".format(jet_idx, postfix),
            expression="jet{}_deepcsv_bb".format(jet_idx),
            binning=binning,
            x_title="Jet_{} prob_{{bb}}".format(jet_idx),
            tags=tags,
        )
        cfg.add_variable(
            name="jet{}_deepcsv_bcomb{}".format(jet_idx, postfix),
            expression="jet{0}_deepcsv_b + jet{0}_deepcsv_bb".format(jet_idx),
            binning=binning,
            x_title="Jet_{} prob_{{b+bb}}".format(jet_idx),
            tags=tags,
        )

# categories
for ch in [ch_ee, ch_emu, ch_mumu]:
    # phase space region loop (measurement, closure, ...)
    for ps_name, ps_sel in get_phasespace_info():
        # inclusive region categories to measure rates
        for rg_name, rg_sel in get_region_info(1, ch, add_btag_cut=False):
            # we skip the emu channel in the lf region because the DY (the main contribution)
            # should have same-flavored leptons
            if rg_name == "LF" and ch == ch_emu:
                continue

            rg_cat_combined = ch.add_category(
                name="{}__{}__{}".format(ch.name, ps_name, rg_name),
                label="{}, {}, {}".format(ch.name, ps_name, rg_name),
                selection=join_root_selection("channel == {}".format(ch.id), ps_sel, rg_sel),
                tags={"scales"},
                aux={
                    "phase_space": ps_name,
                    "region": rg_name,
                }
            )

        # loop over both jet1 jet2 permutations
        for i_tag_jet, i_probe_jet in [(1, 2), (2, 1)]:
            # region loop (hf, lf, ...)
            for rg_name, rg_sel in get_region_info(i_tag_jet, ch):
                if rg_name == "LF" and ch == ch_emu:
                    continue

                rg_cat = ch.add_category(
                    name="{}__{}__{}__j{}".format(ch.name, ps_name, rg_name, i_tag_jet),
                    label="{}, {}, {} region (j{} tagged)".format(ch.name, ps_name, rg_name, i_tag_jet),
                    selection=join_root_selection("channel == {}".format(ch.id), ps_sel, rg_sel),
                )

                # flavor loop (b, c, udsg, ...)
                for fl_name, fl_sel in get_flavor_info(i_probe_jet):
                    fl_cat = rg_cat.add_category(
                        name="{}__f{}".format(rg_cat.name, fl_name),
                        label="{}, {} flavor".format(rg_cat.label, fl_name),
                        selection=join_root_selection(rg_cat.selection, fl_sel),
                    )

                    # pt loop
                    for pt_name, pt_sel, pt_range in get_axis_info(i_probe_jet, "pt")[rg_name]:
                        pt_cat = fl_cat.add_category(
                            name="{}__pt{}".format(fl_cat.name, pt_name),
                            label="{}, pt {}".format(fl_cat.label, pt_name),
                            selection=join_root_selection(fl_cat.selection, pt_sel),
                        )

                        # eta loop
                        for eta_name, eta_sel, eta_range in get_axis_info(i_probe_jet, "abs(eta)", fmt="abs(jet{}_eta)")[rg_name]:
                            eta_cat = pt_cat.add_category(
                                name="{}__eta{}".format(pt_cat.name, eta_name),
                                label="{}, eta {}".format(pt_cat.label, eta_name),
                                selection=join_root_selection(pt_cat.selection, eta_sel),
                                aux={
                                    "channel": ch,
                                    "i_probe_jet": i_probe_jet,
                                    "phase_space": ps_name,
                                    "region": rg_name,
                                    "flavor": fl_name,
                                },
                            )

                            # merged category for both jets and all flavors
                            merged_vars = (ps_name, rg_name, pt_name, eta_name)
                            merged_name = "{}__{}__pt{}__eta{}".format(*merged_vars)
                            if not cfg.has_category(merged_name):
                                label = "{}, {} region, pt {}, eta {}".format(*merged_vars)
                                merged_cat = cfg.add_category(
                                    name=merged_name,
                                    label=label,
                                    tags={"merged"},
                                    aux={
                                        "phase_space": ps_name,
                                        "region": rg_name,
                                        "eta": eta_range,
                                        "pt": pt_range,
                                    }
                                )
                            else:
                                merged_cat = cfg.get_category(merged_name)
                            merged_cat.add_category(eta_cat)

# luminosities per channel in /pb
cfg.set_aux("lumi", {
    #ch_ee: 4767.315,  # B only
    #ch_emu: 4767.315,  # B only
    #ch_mumu: 4767.315,  # B only
    ch_ee: 41296.082,
    ch_emu: 41296.082,
    ch_mumu: 41296.082,
})

# run ranges
rr = cfg.set_aux("run_ranges", {
    "B": (297046, 299329),
    "C": (299368, 302029),
    "D": (302030, 303434),
    "E": (303824, 304797),
    "F": (305040, 306462),
})

# global tags
cfg.set_aux("global_tag", {
    "data": "94X_dataRun2_ReReco_EOY17_v6",
    "mc": "94X_mc2017_realistic_v13",
})

# lumi, normtag and pileup file
cfg.set_aux("lumi_file", "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/"
    "ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt")
cfg.set_aux("normtag_file", "https://www.dropbox.com/s/luj5m5rhb25auhh/normtag_PHYSICS.json")
cfg.set_aux("pileup_file", "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/"
    "PileUp/pileup_latest.txt")

# triggers
cfg.set_aux("triggers", {
    ch_ee: [
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    ],
    ch_emu: [
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    ],
    ch_mumu: [
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",# only 2017B
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*",
    ],
})

# special triggers per real dataset
cfg.set_aux("data_triggers", {
    cfg.get_dataset("data_B_mumu"): [
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
    ],
})

# MET filters
cfg.set_aux("metFilters", {
    "data": [
        "Flag_goodVertices", "Flag_globalTightHalo2016Filter", "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter", "Flag_eeBadScFilter",
        "Flag_ecalBadCalibFilter",
    ],
    "mc": [
        "Flag_goodVertices", "Flag_globalTightHalo2016Filter", "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter", "Flag_ecalBadCalibFilter",
    ],
})

# JER
cfg.set_aux("jer_version", "Summer16_25nsV1")

# JES
cfg.set_aux("jes_version", {
    "data": [
        rr["B"] + ("Fall17_17Nov2017B_V6_DATA",),
        rr["C"] + ("Fall17_17Nov2017C_V6_DATA",),
        rr["D"] + ("Fall17_17Nov2017D_V6_DATA",),
        rr["E"] + ("Fall17_17Nov2017E_V6_DATA",),
        rr["F"] + ("Fall17_17Nov2017F_V6_DATA",),
    ],
    "mc": [
        (1, int(1e9), "Fall17_17Nov2017_V6_MC"),
    ],
})

cfg.set_aux("jes_levels", {
    "data": ["L1FastJet", "L2Relative", "L3Absolute", "L2L3Residual"],
    "mc": ["L1FastJet", "L2Relative", "L3Absolute"],
})

cfg.set_aux("jes_sources", [
    "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation", "SinglePionECAL",
    "SinglePionHCAL", "FlavorQCD", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
    "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal", "RelativeFSR",
    "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF", "PileUpDataMC", "PileUpPtRef",
    "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF",
])

cfg.set_aux("pileup_mc", [
    3.39597497605e-05, 6.63688402133e-06, 1.39533611284e-05, 3.64963078209e-05, 6.00872171664e-05,
    9.33932578027e-05, 0.000120591524486, 0.000128694546198, 0.000361697233219, 0.000361796847553,
    0.000702474896113, 0.00133766053707, 0.00237817050805, 0.00389825605651, 0.00594546732588,
    0.00856825906255, 0.0116627396044, 0.0148793350787, 0.0179897368379, 0.0208723871946,
    0.0232564170641, 0.0249826433945, 0.0262245860346, 0.0272704617569, 0.0283301107549,
    0.0294006137386, 0.0303026836965, 0.0309692426278, 0.0308818046328, 0.0310566806228,
    0.0309692426278, 0.0310566806228, 0.0310566806228, 0.0310566806228, 0.0307696426944,
    0.0300103336052, 0.0288355370103, 0.0273233309106, 0.0264343533951, 0.0255453758796,
    0.0235877272306, 0.0215627588047, 0.0195825559393, 0.0177296309658, 0.0160560731931,
    0.0146022004183, 0.0134080690078, 0.0129586991411, 0.0125093292745, 0.0124360740539,
    0.0123547104433, 0.0123953922486, 0.0124360740539, 0.0124360740539, 0.0123547104433,
    0.0124360740539, 0.0123387597772, 0.0122414455005, 0.011705203844, 0.0108187105305,
    0.00963985508986, 0.00827210065136, 0.00683770076341, 0.00545237697118, 0.00420456901556,
    0.00367513566191, 0.00314570230825, 0.0022917978982, 0.00163221454973, 0.00114065309494,
    0.000784838366118, 0.000533204105387, 0.000358474034915, 0.000238881117601, 0.0001984254989,
    0.000157969880198, 0.00010375646169, 6.77366175538e-05, 4.39850477645e-05, 2.84298066026e-05,
    1.83041729561e-05, 1.17473542058e-05, 7.51982735129e-06, 6.16160108867e-06, 4.80337482605e-06,
    3.06235473369e-06, 1.94863396999e-06, 1.23726800704e-06, 7.83538083774e-07, 4.94602064224e-07,
    3.10989480331e-07, 1.94628487765e-07, 1.57888581037e-07, 1.2114867431e-07, 7.49518929908e-08,
    4.6060444984e-08, 2.81008884326e-08, 1.70121486128e-08, 1.02159894812e-08,
])

cfg.set_aux("min_bias_xs", sn.Number(69.2, (sn.Number.REL, 0.046)))  # mb

# file merging information (stage -> dataset -> files after merging)
cfg.set_aux("file_merging", {
    "trees": {
        "tt_dl": 49,
        "dy_lep_50ToInf_Ht70To100": 5,
        "dy_lep_50ToInf_Ht100To200": 10,
        "dy_lep_50ToInf_Ht200To400": 9,
        "dy_lep_50ToInf_Ht400To600": 10,
        "dy_lep_50ToInf_Ht600To800": 8,
        "dy_lep_50ToInf_Ht800To1200": 3,
        "dy_lep_50ToInf": 7,
    }
})

def get_file_merging(key, dataset):
    dataset_name = dataset if isinstance(dataset, six.string_types) else dataset.name
    return cfg.get_aux("file_merging")[key].get(dataset_name, 1)

cfg.set_aux("get_file_merging", get_file_merging)

# versions
cfg.set_aux("versions", {
    "WriteTrees": "prod3",
    "MergeTrees": "prod3",
    "MergeMetaData": "prod3",
    "WriteHistograms": "prod9",
    "MergeHistograms": "prod9",
    "MeasureScaleFactors": "prod3",
    "FitScaleFactors": "prod3",
})
