# -*- coding: utf-8 -*-

import scinum as sn

def create_config(base_cfg):
    # setup the config for Moriond 2019 (2018 data)
    from analysis.config.campaign_Moriond19 import campaign as campaign_Moriond19
    from analysis.config.jet_tagging_sf import ch_ee, ch_emu, ch_mumu
    cfg = base_cfg.copy(campaign=campaign_Moriond19)

    # add datasets
    dataset_names = [
        "data_A_ee", "data_B_ee", "data_C_ee", "data_D_ee",
        "data_A_emu", "data_B_emu", "data_C_emu", "data_D_emu",
        "data_A_mumu", "data_B_mumu", "data_C_mumu", "data_D_mumu",
        #"data_B_e", "data_C_e", "data_D_e",
        #"data_B_mu", "data_C_mu", "data_D_mu",
        "tt_dl", "tt_sl",
        "dy_lep_10To50",
        "dy_lep_50ToInf",
        #"st_s_lep",
        #"st_t_t", "st_t_tbar",
        "st_tW_t", "st_tW_tbar",
        "WW", "WZ", "ZZ",
        "W_lep",
        "ttH_bb", "ttH_nonbb",
        "ttWJets", "ttZJets",
    ]

    for dataset_name in dataset_names:
        dataset = campaign_Moriond19.get_dataset(dataset_name)
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
            "loose": 0.1241,
            "medium": 0.4184,
            "tight": 0.7527,
        },
        "deepjet": {
            "loose": 0.0494,
            "medium": 0.2770,
            "tight": 0.7264,
        }
    })

    # luminosities per channel in /pb
    cfg.set_aux("lumi", {
        ch_ee: 59970.0, #59966.1613198,
        ch_emu: 59970.0, #59966.1613198,
        ch_mumu: 59970.0, #59966.1613198,
        #ch_e: 59966.1613198,
        #ch_mu: 59966.1613198,
    })

    # run ranges
    # https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis
    rr = cfg.set_aux("run_ranges", {
        "A": (315252, 316995),
        "B": (316998, 319312),
        "C": (319313, 320393),
        "D": (320394, 325273),
    })

    # global tags
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
    cfg.set_aux("global_tag", {
        "data": "102X_dataRun2_Sep2018Rereco_v1", # 102X_dataRun2_Prompt_v11
        "mc": "102X_upgrade2018_realistic_v12",
    })

    # lumi, normtag and pileup file
    cfg.set_aux("lumi_file", "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/"
        "PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt")
    # https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM
    cfg.set_aux("normtag_file", "/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json") # preliminary
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
        #ch_e: [
        #],
        #ch_mu: [
        #],
    })
    # special triggers per real dataset
    cfg.set_aux("data_triggers", {})

    # MET filters
    # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2?rev=131
    cfg.set_aux("metFilters", {
        "data": [
            "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter",
            "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter", "Flag_eeBadScFilter",
            "Flag_ecalBadCalibReducedMINIAODFilter", # TODO: Implement updated filter
        ],
        "mc": [
            "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter",
            "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter",
            "Flag_ecalBadCalibReducedMINIAODFilter",
        ],
    })

    # JER
    # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    cfg.set_aux("jer_version", "Fall17_V3")

    # JES
    # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
    cfg.set_aux("jes_version", {
        "data": [
            rr["A"] + ("Fall17_09May2018F_V3_DATA",),
            rr["B"] + ("Fall17_09May2018F_V3_DATA",),
            rr["C"] + ("Fall17_09May2018F_V3_DATA",),
            rr["D"] + ("Fall17_09May2018F_V3_DATA",),
        ],
        "mc": [
            (1, int(1e9), "Fall17_17Nov2017_V32_MC"),
        ],
    })

    # https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2018_25ns_JuneProjectionFull18_PoissonOOTPU_cfi.py
    cfg.set_aux("pileup_mc", [
        4.695341e-10, 1.206213e-06, 1.162593e-06, 6.118058e-06, 1.626767e-05,
        3.508135e-05, 7.12608e-05, 0.0001400641, 0.0002663403, 0.0004867473,
        0.0008469, 0.001394142, 0.002169081, 0.003198514, 0.004491138,
        0.006036423, 0.007806509, 0.00976048, 0.0118498, 0.01402411,
        0.01623639, 0.01844593, 0.02061956, 0.02273221, 0.02476554,
        0.02670494, 0.02853662, 0.03024538, 0.03181323, 0.03321895,
        0.03443884, 0.035448, 0.03622242, 0.03674106, 0.0369877,
        0.03695224, 0.03663157, 0.03602986, 0.03515857, 0.03403612,
        0.0326868, 0.03113936, 0.02942582, 0.02757999, 0.02563551,
        0.02362497, 0.02158003, 0.01953143, 0.01750863, 0.01553934,
        0.01364905, 0.01186035, 0.01019246, 0.008660705, 0.007275915,
        0.006043917, 0.004965276, 0.004035611, 0.003246373, 0.002585932,
        0.002040746, 0.001596402, 0.001238498, 0.0009533139, 0.0007282885,
        0.000552306, 0.0004158005, 0.0003107302, 0.0002304612, 0.0001696012,
        0.0001238161, 8.96531e-05, 6.438087e-05, 4.585302e-05, 3.23949e-05,
        2.271048e-05, 1.580622e-05, 1.09286e-05, 7.512748e-06, 5.140304e-06,
        3.505254e-06, 2.386437e-06, 1.625859e-06, 1.111865e-06, 7.663272e-07,
        5.350694e-07, 3.808318e-07, 2.781785e-07, 2.098661e-07, 1.642811e-07,
        1.312835e-07, 1.081326e-07, 9.141993e-08, 7.890983e-08, 6.91468e-08,
        6.119019e-08, 5.443693e-08, 4.85036e-08, 4.31486e-08, 3.822112e-08
    ])

    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II
    cfg.set_aux("min_bias_xs", sn.Number(69.2, (sn.Number.REL, 0.046)))  # mb

    # file merging information (stage -> dataset -> files after merging)
    cfg.set_aux("file_merging", {
        "trees": {
            "tt_dl": 88,
            "tt_sl": 2,
            "dy_lep_50ToInf": 19,
            "st_tW_t": 2,
            "ttH_bb": 3,
            "ttH_nonbb": 3,
            "ttWJets": 5,
            "ttZJets": 11,
        }
    })

    # versions
    cfg.set_aux("versions", {
        "WriteTrees": "prod2",
        "MergeTrees": "prod2",
        "MergeMetaData": "prod1",
        "WriteHistograms": "test2",
        "MergeHistograms": "test2",
        "MeasureCScaleFactors": "test1",
        "MeasureScaleFactors": "test1",
        "FitScaleFactors": "test1",
        "GetScaleFactorWeights": "test1",
        "MergeScaleFactorWeights": "test1",
    })

    return cfg
