# -*- coding: utf-8 -*-

import scinum as sn
import numpy as np

def create_config(base_cfg):
    # setup the config for 2018 data
    from analysis.config.campaign_UltraLegacy18 import campaign as campaign_UltraLegacy18
    from analysis.config.jet_tagging_sf import ch_ee, ch_emu, ch_mumu, ch_e, ch_mu
    cfg = base_cfg.copy(campaign=campaign_UltraLegacy18)

    # add datasets
    dataset_names = [
        "data_A_ee", "data_B_ee", "data_C_ee", "data_D_ee",
        "data_A_emu", "data_B_emu", "data_C_emu", "data_D_emu",
        "data_A_mumu", "data_B_mumu", "data_C_mumu", "data_D_mumu",
        "data_A_e", "data_B_e", "data_C_e", "data_D_e",
        "data_A_mu", "data_B_mu", "data_C_mu", "data_D_mu",
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
        dataset = campaign_UltraLegacy18.get_dataset(dataset_name)
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
            "loose": 0.1208,
            "medium": 0.4168,
            "tight": 0.7665,
        },
        "deepjet": {
            "loose": 0.0490,
            "medium": 0.2783,
            "tight": 0.7100,
        }
    })

    # luminosities per channel in /pb
    cfg.set_aux("lumi", {
        ch_ee: 59830.,
        ch_emu: 59830.,
        ch_mumu: 59830.,
        ch_e: 59830.,
        ch_mu: 59830.,
    })

    # run ranges
    rr = cfg.set_aux("run_ranges", {
        "A": (315252, 316995),
        "B": (316998, 319312),
        "C": (319313, 320393),
        "D": (320394, 325273),
    })

    # global tags
    cfg.set_aux("global_tag", {
        "data": "106X_dataRun2_v28",
        "mc": "106X_upgrade2018_realistic_v11_L1v1",
    })

    # lumi, normtag and pileup file
    cfg.set_aux("lumi_file", "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/"
        "Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt")
    # https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM
    cfg.set_aux("normtag_file", "/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json")
    cfg.set_aux("pileup_file", "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/"
        "PileUp/pileup_latest.txt")

    # triggers
    # https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2018
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
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*",
        ],
	ch_e: [
            "HLT_Ele35_WPTight_Gsf_v*",
            "HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v*",
        ],
	ch_mu: [
            "HLT_IsoMu24_v*",
        ],
    })
    # special triggers per real dataset
    cfg.set_aux("data_triggers", {})


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
    cfg.set_aux("jer_version", "Summer19UL18_JRV2")

    # JES
    cfg.set_aux("jes_version", {
        "data": [
            rr["A"] + ("Summer19UL18_RunA_V5_DATA",),
            rr["B"] + ("Summer19UL18_RunB_V5_DATA",),
            rr["C"] + ("Summer19UL18_RunC_V5_DATA",),
            rr["D"] + ("Summer19UL18_RunD_V5_DATA",),
        ],
        "mc": [
            (1, int(1e9), "Summer19UL18_V5_MC"),
        ],
    })
    # JES veto maps
    cfg.set_aux("jes_veto_map", {
        "file": "Summer19UL18_V1/hotjets-UL18.root",
        "hist_name": "h2hot_ul18_plus_hem1516_plus_hbp2m1",
    })

    cfg.set_aux("jes_uncertainty_file", {
        "factorized": None,  # take file from jes github
        "reduced": "",
    })

    # https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2018_25ns_UltraLegacy_PoissonOOTPU_cfi.py
    cfg.set_aux("pileup_mc", [
        8.89374611122e-07, 1.1777062868e-05, 3.99725585118e-05, 0.000129888015252, 0.000265224848687,
        0.000313088635109, 0.000353781668514, 0.000508787237162, 0.000873670065767, 0.00147166880932,
        0.00228230649018, 0.00330375581273, 0.00466047608406, 0.00624959203029, 0.00810375867901,
        0.010306521821, 0.0129512453978, 0.0160303925502, 0.0192913204592, 0.0223108613632,
        0.0249798930986, 0.0273973789867, 0.0294402350483, 0.031029854302, 0.0324583524255,
        0.0338264469857, 0.0351267479019, 0.0360320204259, 0.0367489568401, 0.0374133183052,
        0.0380352633799, 0.0386200967002, 0.039124376968, 0.0394201612616, 0.0394673457109,
        0.0391705388069, 0.0384758587461, 0.0372984548399, 0.0356497876549, 0.0334655175178,
        0.030823567063, 0.0278340752408, 0.0246009685048, 0.0212676009273, 0.0180250593982,
        0.0149129830776, 0.0120582333486, 0.00953400069415, 0.00738546929512, 0.00563442079939,
        0.00422052915668, 0.00312446316347, 0.00228717533955, 0.00164064894334, 0.00118425084792,
        0.000847785826565, 0.000603466454784, 0.000419347268964, 0.000291768785963, 0.000199761337863,
        0.000136624574661, 9.46855200945e-05, 6.80243180179e-05, 4.94806013765e-05, 3.53122628249e-05,
        2.556765786e-05, 1.75845711623e-05, 1.23828210848e-05, 9.31669724108e-06, 6.0713272037e-06,
        3.95387384933e-06, 2.02760874107e-06, 1.22535149516e-06, 9.79612472109e-07, 7.61730246474e-07,
        4.2748847738e-07, 2.41170461205e-07, 1.38701083552e-07, 3.37678010922e-08, 0.0,
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
            "data_D_e": 2,
            "data_A_mu": 2,
            "data_D_mu": 3,
            "tt_dl": 456,
            "tt_sl": 491,
            "dy_lep_LO_50ToInf": 30,
            "st_s_lep": 14,
            "st_t_t": 14,
            "st_t_tbar": 7,
            "st_tW_t": 34,
            "st_tW_tbar": 31,
            "WW": 3,
            "WZ": 2,
            "W_lep": 3
        }
    })

    # versions
    cfg.set_aux("versions", {
        "WriteTrees": "prod2", # including SL events
        "MergeTrees": "prod2",
        "MergeMetaData": "prod2",
        "WriteHistograms": "prod2",
        "MergeHistograms": "prod2",
        "MeasureCScaleFactors": "prod1",
        "MeasureScaleFactors": "prod1",
        "FitScaleFactors": "prod1",
        "BundleScaleFactors": "prod1",
        "GetScaleFactorWeights": "prod1",
        "MergeScaleFactorWeights": "prod1",
        "OptimizeBinning": "prod1",
        "CreateScaleFactorResults": "prod1",
    })

    return cfg
