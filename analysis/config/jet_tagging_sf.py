# -*- coding: utf-8 -*-
# flake8: noqa

"""
Definition of the analysis for extracting jet tagging scale factors
as well as its configurations for different campaigns.
"""


import re

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

# store DeepCSV working points

cfg.set_aux("working_points", {
    "deepcsv": {
        "loose": 0.1522,
        "medium": 0.4941,
        "tight": 0.8001,
    }
})

# link processes
cfg.add_process(process_data_ee)
cfg.add_process(process_data_emu)
cfg.add_process(process_data_mumu)
cfg.add_process(process_tt_dl)
# cfg.add_process(process_dy_lep)  # tmp
cfg.add_process(process_st_tW)
# cfg.add_process(process_WW_sl)  # tmp

# add datasets
# dataset_names = [  # tmp
#     "data_B_ee", "data_C_ee", "data_D_ee", "data_E_ee", "data_F_ee",
#     "data_B_emu", "data_C_emu", "data_D_emu", "data_E_emu", "data_F_emu",
#     "data_B_mumu", "data_C_mumu", "data_D_mumu", "data_E_mumu", "data_F_mumu",
#     "tt_dl", "dy_lep_5To50_Ht70To100", "dy_lep_5To50_Ht100To200", "dy_lep_5To50_Ht200To400",
#     "dy_lep_5To50_Ht600ToInf", "dy_lep_50ToInf_Ht70To100", "dy_lep_50ToInf_Ht100To200",
#     "dy_lep_50ToInf_Ht200To400", "dy_lep_50ToInf_Ht400To600", "dy_lep_50ToInf_Ht600To800",
#     "dy_lep_50ToInf_Ht800To1200", "dy_lep_50ToInf_Ht1200To2500", "dy_lep_50ToInf_Ht2500ToInf",
#     "st_tW_t", "st_tW_tbar", "WW_sl",
# ]
dataset_names = [
    "data_B_ee", "data_B_emu", "data_B_mumu", "tt_dl", "st_tW_t",
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

# define nested categories (analysis phase space -> hf/lf region -> flavor -> pt bin -> eta bin)
def get_phasespace_info():
    return [
        ("measure", join_root_selection("n_jets == 2")),
        ("closure", join_root_selection("n_jets >= 2")),
    ]

def get_region_info(idx, channel, btagger="deepcsv_bcomb", et_miss=30., z_window=10.):
    csv_tight = cfg.get_aux("working_points")["deepcsv"]["tight"]
    csv_loose = cfg.get_aux("working_points")["deepcsv"]["loose"]

    hf_cuts, lf_cuts = [], []
    # jet tagging requirement
    hf_cuts.append("jet{}_{} > {}".format(btagger, idx, csv_tight))
    lf_cuts.append("jet{}_{} < {}".format(btagger, idx, csv_loose))

    # the following cuts do not apply for emu
    if channel != "emu":
        # ET-miss requirement
        hf_cuts.append(" > {}".format(et_miss))
        lf_cuts.append(" < {}".format(et_miss))

        # z-mass window
        hf_cuts.append("abs(mll - {}) > {}".format(Z_MASS.nominal, z_window))
        lf_cuts.append("abs(mll - {}) < {}".format(Z_MASS.nominal, z_window))

        # z peak diamond
        # lf_cuts.append("pass_z_peak == 0")

    return [
        ("HF", join_root_selection(hf_cuts)),
        ("LF", join_root_selection(lf_cuts)),
    ]

def get_flavor_info(idx):
    return [
        ("b", "abs(jet%d_flavor) == 5" % idx),
        ("c", "abs(jet%d_flavor) == 4" % idx),
        ("udsg", "abs(jet%d_flavor) != 5 && abs(jet%d_flavor) != 4" % (idx, idx)),
    ]

def get_pt_info(idx):
    return {
        "b": [
            ("20To30", "jet%d_pt > 20 && jet%d_pt <= 30" % (idx, idx)),
            ("30To50", "jet%d_pt > 30 && jet%d_pt <= 50" % (idx, idx)),
            ("50To70", "jet%d_pt > 50 && jet%d_pt <= 70" % (idx, idx)),
            ("70To100", "jet%d_pt > 70 && jet%d_pt <= 100" % (idx, idx)),
            ("100ToInf", "jet%d_pt > 100" % idx),
        ],
        "c": [
            ("20To30", "jet%d_pt > 20 && jet%d_pt <= 30" % (idx, idx)),
            ("30To50", "jet%d_pt > 30 && jet%d_pt <= 50" % (idx, idx)),
            ("50To70", "jet%d_pt > 50 && jet%d_pt <= 70" % (idx, idx)),
            ("70To100", "jet%d_pt > 70 && jet%d_pt <= 100" % (idx, idx)),
            ("100ToInf", "jet%d_pt > 100" % idx),
        ],
        "udsg": [
            ("20To30", "jet%d_pt > 20 && jet%d_pt <= 30" % (idx, idx)),
            ("30To40", "jet%d_pt > 30 && jet%d_pt <= 40" % (idx, idx)),
            ("40To60", "jet%d_pt > 40 && jet%d_pt <= 60" % (idx, idx)),
            ("60ToInf", "jet%d_pt > 60" % idx),
        ]
    }

def get_eta_info(idx):
    return {
        "b": [
            ("0To2p4", "jet%d_eta <= 2.4" % idx),
        ],
        "c": [
            ("0To2p4", "jet%d_eta <= 2.4" % idx),
        ],
        "udsg": [
            ("0To0p8", "jet%d_eta <= 0.8" % idx),
            ("0p8To1p6", "jet%d_eta > 0.8 && jet%d_eta <= 1.6" % (idx, idx)),
            ("1p6To2p4", "jet%d_eta > 1.6 && jet%d_eta <= 2.4" % (idx, idx)),
        ]
    }

for ch in [ch_ee, ch_emu, ch_mumu]:
    # phase space region loop (measurement, closure, ...)
    for ps_name, ps_sel in get_phasespace_info():
        # loop over both jet1 jet2 permutations
        for i_tag_jet, i_probe_jet in [(1, 2), (2, 1)]:
            # region loop (hf, lf, ...)
            for rg_name, rg_sel in get_region_info(i_tag_jet, ch):
                rg_cat = ch.add_category(
                    name="{}__{}__{}__j{}".format(ch.name, ps_name, rg_name, i_tag_jet),
                    label="{}, {} region (j{} tagged)".format(ch.name, ps_name, i_tag_jet),
                    selection=join_root_selection("channel == {}".format(ch.id), ps_sel, rg_sel),
                )
                # flavor loop (b, c, udsg, ...)
                for fl_name, fl_sel in get_flavor_info(i_probe_jet):
                    fl_cat = rg_cat.add_category(
                        name="{}__{}".format(rg_cat.name, fl_name),
                        label="{}, {} flavor".format(rg_cat.label, fl_name),
                        selection=join_root_selection(rg_cat.selection, fl_sel),
                    )
                    # pt loop
                    for pt_name, pt_sel in get_pt_info(i_probe_jet)[fl_name]:
                        pt_cat = fl_cat.add_category(
                            name="{}__pt{}".format(fl_cat.name, pt_name),
                            label="{}, pt {}".format(fl_cat.label, pt_name),
                            selection="{} && {}".format(fl_cat.selection, pt_sel),
                        )
                        # eta loop
                        for eta_name, eta_sel in get_eta_info(i_probe_jet)[fl_name]:
                            eta_cat = pt_cat.add_category(
                                name="{}__eta{}".format(pt_cat.name, eta_name),
                                label="{}, eta {}".format(pt_cat.label, eta_name),
                                selection="{} && {}".format(pt_cat.selection, eta_sel),
                            )

                            # merged category for both jets
                            merged_name = re.sub(r"__j\d+__", "__", eta_cat.name)
                            if not ch.has_category(merged_name):
                                merged_cat = ch.add_category(
                                    name=merged_name,
                                    label=re.sub(r" \(j\d+ tagged\)", "", eta_cat.label),
                                )
                            else:
                                merged_cat = ch.get_category(merged_name)
                            merged_cat.add_category(eta_cat)

# variables
cfg.add_variable(
    name="jet1_pt",
    expression="jet1_pt",
    binning=(25, 0., 500.,),
    unit="GeV",
    x_title="Jet_{1} p_{T}",
)
cfg.add_variable(
    name="jet1_deepcsv_b",
    expression="jet1_deepcsv_b",
    binning=(25, 0., 1.,),
    x_title="Jet_{1} prob_{b}",
)
cfg.add_variable(
    name="jet1_deepcsv_bb",
    expression="jet1_deepcsv_bb",
    binning=(25, 0., 1.,),
    x_title="Jet_{1} prob_{bb}",
)
cfg.add_variable(
    name="jet1_deepcsv_bcomb",
    expression="jet1_deepcsv_b + jet1_deepcsv_bb",
    binning=(25, 0., 1.,),
    x_title="Jet_{1} prob_{b+bb}",
)

# luminosities per channel in /pb
cfg.set_aux("lumi", {
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

# lumi and normtag file
cfg.set_aux("lumi_file", "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/"
    "ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt")
cfg.set_aux("normtag_file", "https://www.dropbox.com/s/luj5m5rhb25auhh/normtag_PHYSICS.json")

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
        # note: the 2017B mumu dataset uses a different trigger
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

# file merging information (stage -> dataset -> files after merging)
cfg.set_aux("file_merging", {
    "trees": {
        "tt_dl": 74,
    }
})

def get_file_merging(key, dataset):
    dataset_name = dataset if isinstance(dataset, six.string_types) else dataset.name
    return cfg.get_aux("file_merging")[key].get(dataset_name, 1)

cfg.set_aux("get_file_merging", get_file_merging)

# versions
cfg.set_aux("versions", {
    "WriteTrees": "prod1",
    "MergeTrees": "prod1",
    "MergeMetaData": "prod1",
    "WriteHistograms": "prod1",
})
