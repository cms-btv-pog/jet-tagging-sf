# -*- coding: utf-8 -*-


__all__ = ["campaign", "analysis"]


import order as od
import scinum as sn


# constants
BR_W_HAD = sn.Number(0.6741, {"br_w": 0.0027})
BR_W_LEP = 1 - BR_W_HAD
BR_WW_SL = 2 * BR_W_HAD.mul(BR_W_LEP, rho=-1, inplace=False)
BR_WW_DL = BR_W_LEP**2
BR_WW_FH = BR_W_HAD**2

# campaign
campaign = od.Campaign("2017_13Tev_25ns", 1, ecm=13000, bx=25)

# processes
process_data_ee = od.Process(
    "data_ee", 1,
    is_data=True,
    label="data",
)

process_tt = od.Process(
    "tt", 10,
    label=r"$t\bar{t}$ + Jets",
    xsecs={
        13: sn.Number(831.76, {
            "scale": (19.77, 29.20),
            "pdf": 35.06,
            "mtop": (23.18, 22.45),
        }),
    },
)

process_tt_dl = od.Process(
    "tt_dl", 12,
    label=r"$t\bar{t}$ + Jets, DL",
    xsecs={
        13: process_tt.get_xsec(13) * BR_WW_DL,
    },
)

# datasets
dataset_data_ee = od.Dataset(
    "data_ee", 1,
    campaign=campaign,
    is_data=True,
    n_files=759,
    n_events=58088760,
    keys=["/DoubleEG/Run2017B-17Nov2017-v1/MINIAOD"],
)

dataset_tt_dl = od.Dataset(
    "tt_dl", 12,
    campaign=campaign,
    n_files=164 + 1230,
    n_events=8705576 + 69705626,
    keys=[
        "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/MINIAODSIM",
        "/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM",
    ],
)

# link processes to datasets
dataset_data_ee.add_process(process_data_ee)
dataset_tt_dl.add_process(process_tt_dl)

# add the analysis and a config for the 2017 campaign
analysis = od.Analysis("jet-tagging-sf", 1)
cfg = analysis.add_config(campaign=campaign)

# link processes
cfg.add_process(process_data_ee)
cfg.add_process(process_tt_dl)

# add datasets
cfg.add_dataset(dataset_data_ee)
cfg.add_dataset(dataset_tt_dl)

# define channels
ch_ee = cfg.add_channel("ee", 1)
ch_emu = cfg.add_channel("emu", 2)
ch_mumu = cfg.add_channel("mumu", 3)

# store channels per real dataset
cfg.set_aux("dataset_channels", {
    dataset_data_ee: ch_ee,
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
cfg.set_aux("lumi_file", "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt")
cfg.set_aux("normtag_file", "/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json")

# triggers
cfg.set_aux("triggers", {
    ch_ee: [
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    ],
    ch_emu: [
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
    ],
    ch_mumu: [
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
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
cfg.set_aux("jer_version", "Spring16_25nsV10")

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
