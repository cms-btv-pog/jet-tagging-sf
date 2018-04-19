# -*- coding: utf-8 -*-


__all__ = ["campaign", "analysis"]


import order as od
import scinum as sn


# campaign
campaign = od.Campaign("2017_13Tev_25ns", 1, ecm=13000, bx=25)

# processes
process_data_ee = od.Process(
    "data_ee", 1,
    is_data=True,
    label="data",
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

# link processes to datasets
dataset_data_ee.add_process(process_data_ee)


# add the analysis and a config for the 2017 campaign
analysis = od.Analysis("jet-tagging-sf", 1)
cfg = analysis.add_config(campaign=campaign)

# link processes
cfg.add_process(process_data_ee)

# add datasets
cfg.add_dataset(dataset_data_ee)

# define channels
ch_ee = cfg.add_channel("ee", 1)
ch_emu = cfg.add_channel("emu", 2)
ch_mumu = cfg.add_channel("mumu", 3)

# store channels per real dataset
cfg.set_aux("dataset_channels", {
    dataset_data_ee: ch_ee,
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
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    ],
    ch_mumu: [
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
    ],
})
