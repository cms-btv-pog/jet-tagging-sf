# -*- coding: utf-8 -*-

"""
Definition of the analysis for extracting jet tagging scale factors
as well as its configurations for different campaigns.
"""


import order as od
import scinum as sn

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
# cfg.add_process(process_dy_lep)  # tmp
cfg.add_process(process_st_tW)
# cfg.add_process(process_WW_sl)  # tmp

# add datasets
# dataset_names = [  # tmp
#     "data_B_ee", "data_C_ee", "data_D_ee", "data_E_ee", "data_F_ee",
#     "data_B_emu", "data_C_emu", "data_D_emu", "data_E_emu", "data_F_emu",
#     "data_B_mumu", "data_C_mumu", "data_D_mumu", "data_E_mumu", "data_F_mumu",
#     "tt_dl", "dy_lep_10To50", "dy_lep_50ToInf", "st_tW_t", "st_tW_tbar", "WW_sl",
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
# TODO: this should be measured, see #5
cfg.set_aux("file_merging", {
    "trees": {
        "data_B_ee": 1,
        "data_B_emu": 1,
        "data_B_mumu": 1,
        "st_tW_t": 1,
    }
})


# versions
cfg.set_aux("versions", {
    "WriteTrees": "prod1",
})
