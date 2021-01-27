# -*- coding: utf-8 -*-

import scinum as sn
import numpy as np

def create_config(base_cfg):
    # setup the config for 2017 data
    from analysis.config.campaign_UltraLegacy17 import campaign as campaign_UltraLegacy17
    from analysis.config.jet_tagging_sf import ch_ee, ch_emu, ch_mumu, ch_e, ch_mu
    cfg = base_cfg.copy(campaign=campaign_UltraLegacy17)

    # add datasets
    dataset_names = [
        "data_B_ee", "data_C_ee", "data_D_ee", "data_E_ee", "data_F_ee",
        "data_B_emu", "data_C_emu", "data_D_emu", "data_E_emu", "data_F_emu",
        "data_B_mumu", "data_C_mumu", "data_D_mumu", "data_E_mumu", "data_F_mumu",
        "data_B_e", "data_C_e", "data_D_e", "data_E_e", "data_F_e",
        "data_B_mu", "data_C_mu", "data_D_mu", "data_E_mu", "data_F_mu",
        "tt_dl", "tt_sl",
        "dy_lep_10To50",
        #"dy_lep_50ToInf",
        "dy_lep_LO_50ToInf",
        #"dy_lep_0Jets", "dy_lep_1Jets", "dy_lep_2Jets",
        "st_s_lep",
        "st_t_t", "st_t_tbar",
        "st_tW_t", "st_tW_tbar",
        "WW", "WZ", "ZZ",
        "W_lep",
        #"ttH",
        #"ttWJets_lep", "ttWJets_had", "ttZJets_lep", "ttZJets_had",
    ]

    for dataset_name in dataset_names:
        dataset = campaign_UltraLegacy17.get_dataset(dataset_name)
        cfg.add_dataset(dataset)

    # store channels per real dataset
    cfg.set_aux("dataset_channels", {
        dataset: cfg.get_channel(dataset.name.split("_")[-1])
        for dataset in cfg.datasets.values()
        if dataset.is_data
    })

    # store b-tagger working points
    cfg.set_aux("working_points", {
        "deepcsv": {
            "loose": 0.1355,
            "medium": 0.4506,
            "tight": 0.7738,
        },
        "deepjet": {
            "loose": 0.0532,
            "medium": 0.3040,
            "tight": 0.7476,
        }
    })

    # luminosities per channel in /pb
    cfg.set_aux("lumi", {
        ch_ee: 41480.0,
        ch_emu: 41480.0,
        ch_mumu: 41480.0,
        ch_e: 41480.0,
        ch_mu: 41480.0,
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
        "data": "106X_dataRun2_v28",
        "mc": "106X_mc2017_realistic_v7",
    })

    # lumi, normtag and pileup file
    cfg.set_aux("lumi_file", "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/"
        "13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt")
    # https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM
    cfg.set_aux("normtag_file", "/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json")
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
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*", # only 2017 C-F
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*", # only 2017 C-F
        ],
        ch_e: [
            "HLT_Ele35_WPTight_Gsf_v*",
            "HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v*",
        ],
        ch_mu: [
            "HLT_IsoMu27_v*",
            "HLT_IsoMu24_eta2p1_v*", # only 2017 B, C, D
        ],
    })

    # special triggers per real dataset
    #cfg.set_aux("data_triggers", {
    #    cfg.get_dataset("data_B_mumu"): [
    #        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
    #    ],
    #})
    #for era in ["C", "D", "E", "F"]:
    #    cfg.set_aux("data_triggers", {
    #        cfg.get_dataset("data_{}_mumu".format(era)): [
    #            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",
    #            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*",
    #        ],
    #    })
    for era in ["E", "F"]:
        cfg.set_aux("data_triggers", {
            cfg.get_dataset("data_{}_mu".format(era)): [
                "HLT_IsoMu27_v*",
            ],
        })

    # MET filters
    # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    cfg.set_aux("metFilters", {
        "data": [
            "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter",
            "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter", #"Flag_BadChargedCandidateFilter",
            "Flag_eeBadScFilter", #"Flag_ecalBadCalibReducedMINIAODFilter",
        ],
        "mc": [
            "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter",
            "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter", #"Flag_BadChargedCandidateFilter",
            #"Flag_ecalBadCalibReducedMINIAODFilter",
        ],
    })

    # JER
    cfg.set_aux("jer_version", "Summer19UL17")

    # JES
    cfg.set_aux("jes_version", {
        "data": [
            rr["B"] + ("Summer19UL17_RunB_V5_DATA",),
            rr["C"] + ("Summer19UL17_RunC_V5_DATA",),
            rr["D"] + ("Summer19UL17_RunD_V5_DATA",),
            rr["E"] + ("Summer19UL17_RunE_V5_DATA",),
            rr["F"] + ("Summer19UL17_RunF_V5_DATA",),
        ],
        "mc": [
            (1, int(1e9), "Summer19UL17_V5_MC"),
        ],
    })

    cfg.set_aux("jes_uncertainty_file", {
        "factorized": None,  # take file from jes github
        "reduced": "https://cernbox.cern.ch/index.php/s/4ks9YSbEc7ietyE/download",
    })

    # https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2017_25ns_UltraLegacy_PoissonOOTPU_cfi.py
    cfg.set_aux("pileup_mc", [
    1.1840841518e-05, 3.46661037703e-05, 8.98772521472e-05, 7.47400487733e-05, 0.000123005176624,
    0.000156501700614, 0.000154660478659, 0.000177496185603, 0.000324149805611, 0.000737524009713,
    0.00140432980253, 0.00244424508696, 0.00380027898037, 0.00541093042612, 0.00768803501793,
    0.010828224552, 0.0146608623707, 0.01887739113, 0.0228418813823, 0.0264817796874,
    0.0294637401336, 0.0317960986171, 0.0336645950831, 0.0352638818387, 0.036869429333,
    0.0382797316998, 0.039386705577, 0.0398389681346, 0.039646211131, 0.0388392805703,
    0.0374195678161, 0.0355377892706, 0.0333383902828, 0.0308286549265, 0.0282914440969,
    0.0257860718304, 0.02341635055, 0.0213126338243, 0.0195035612803, 0.0181079838989,
    0.0171991315458, 0.0166377598339, 0.0166445341361, 0.0171943735369, 0.0181980997278,
    0.0191339792146, 0.0198518804356, 0.0199714909193, 0.0194616474094, 0.0178626975229,
    0.0153296785464, 0.0126789254325, 0.0100766041988, 0.00773867100481, 0.00592386091874,
    0.00434706240169, 0.00310217013427, 0.00213213401899, 0.0013996000761, 0.000879148859271,
    0.000540866009427, 0.000326115560156, 0.000193965828516, 0.000114607606623, 6.74262828734e-05,
    3.97805301078e-05, 2.19948704638e-05, 9.72007976207e-06, 4.26179259146e-06, 2.80015581327e-06,
    1.14675436465e-06, 2.52452411995e-07, 9.08394910044e-08, 1.14291987912e-08, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0
    ])

    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II
    cfg.set_aux("min_bias_xs", sn.Number(69.2, (sn.Number.REL, 0.046)))  # mb

    # file merging information (stage -> dataset -> files after merging)
    cfg.set_aux("file_merging", {
        "trees": {
#            "tt_dl": 155,
#            "tt_sl": 4,
#            "dy_lep_50ToInf": 73,
#            "dy_lep_LO_50ToInf": 21,
#            "st_tW_t": 4,
#            "st_tW_tbar": 4,
            "tt_dl": 194,
            "tt_sl": 540,
            "dy_lep_LO_50ToInf": 8,
            "st_s_lep": 14,
            "st_t_t": 2,
            "st_tW_t": 32,
            "st_tW_tbar": 25,
            "WW": 2,
            "W_lep": 4,
        }
    })

    # versions
    cfg.set_aux("versions", {
        "WriteTrees": "prod2_sl",
        "MergeTrees": "prod2_sl",
        "MergeMetaData": "prod2_sl",
        "WriteHistograms": "prod5_sl",
        "MergeHistograms": "prod5_sl",
        "MeasureCScaleFactors": "prod5",
        "MeasureScaleFactors": "prod5",
        "FitScaleFactors": "prod5",
        "BundleScaleFactors": "prod5",
        "GetScaleFactorWeights": "prod5",
        "MergeScaleFactorWeights": "prod5",
        "OptimizeBinning": "prod4",
        "CreateScaleFactorResults": "prod5",
    })

    return cfg
