# -*- coding: utf-8 -*-

import scinum as sn

def create_config(base_cfg):
    # setup the config for Moriond 2019 legacy (2016 data)
    from analysis.config.campaign_Moriond19_legacy import campaign as campaign_Moriond19_legacy
    from analysis.config.jet_tagging_sf import ch_ee, ch_emu, ch_mumu
    cfg = base_cfg.copy(campaign=campaign_Moriond19_legacy)

    # add datasets
    dataset_names = [
        "data_B_ee", "data_C_ee", "data_D_ee", "data_E_ee", "data_F_ee",  "data_G_ee", "data_H_ee",
        "data_B_emu", "data_C_emu", "data_D_emu", "data_E_emu", "data_F_emu", "data_G_emu", "data_H_emu",
        "data_B_mumu", "data_C_mumu", "data_D_mumu", "data_E_mumu", "data_F_mumu", "data_G_mumu", "data_H_mumu",
        #"data_B_e", "data_C_e", "data_D_e", "data_E_e", "data_F_e", "data_G_e", "data_H_e",
        #"data_B_mu", "data_C_mu", "data_D_mu", "data_E_mu", "data_F_mu", "data_G_mu", "data_H_mu",
        "tt_dl", "tt_sl",
        "dy_lep_10To50",
        "dy_lep_50ToInf",
        "st_s_lep",
        "st_t_t", "st_t_tbar",
        "st_tW_t", "st_tW_tbar",
        "WW", "WZ", "ZZ",
        "W_lep",
        #"ttH_bb", "ttH_nonbb",
        "ttWJets_lep", "ttWJets_had", "ttZJets_lep", "ttZJets_had",
    ]

    for dataset_name in dataset_names:
        dataset = campaign_Moriond19_legacy.get_dataset(dataset_name)
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
            "loose": 0.2217,
            "medium": 0.6321,
            "tight": 0.8953,
        },
        "deepflavor": {
            "loose": 0.0614,
            "medium": 0.3093,
            "tight": 0.7221,
        }
    })

    # luminosities per channel in /pb
    cfg.set_aux("lumi", { # TODO: Calculate exact number in task
        ch_ee: 35922.0,
        ch_emu: 35922.0,
        ch_mumu: 35922.0,
        #ch_e: 35922.0,
        #ch_mu: 35922.0,
    })

    # run ranges
    # https://twiki.cern.ch/CMS/PdmV2016Analysis#Datasets_to_be_used
    rr = cfg.set_aux("run_ranges", {
        "B": (272007, 275376),
        "C": (275657, 276283),
        "D": (276315, 276811),
        "E": (276831, 277420),
        "F": (277772, 278808),
        "G": (278820, 280385),
        "H": (280919, 284044),
    })

    # global tags
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
    cfg.set_aux("global_tag", {
        "data": "94X_dataRun2_v10",
        "mc": "94X_mcRun2_asymptotic_v3",
    })

    # lumi, normtag and pileup file
    cfg.set_aux("lumi_file", "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/"
        "Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt")
    # https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM
    cfg.set_aux("normtag_file", "/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json")
    cfg.set_aux("pileup_file", "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/"
        "PileUp/pileup_latest.txt")

    # triggers
    cfg.set_aux("triggers", {
        ch_ee: [
            "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", # TODO: Reduced json?
        ],
        ch_emu: [
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        ],
        ch_mumu: [
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
        ],
        #ch_e: [
        #],
        #ch_mu: [
        #],
    })
    # special triggers per real dataset
    cfg.set_aux("data_triggers", {})

    # MET filters
    cfg.set_aux("metFilters", {
    # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2?rev=131
        "data": [
            "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter",
            "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter", "Flag_eeBadScFilter",
        ],
        "mc": [
            "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter",
            "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter",
        ],
    })

    # JER
    # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    cfg.set_aux("jer_version", "Summer16_25nsV1")

    # JES
    # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
    cfg.set_aux("jes_version", {
        "data": [
            rr["B"] + ("Summer16_07Aug2017BCD_V11_DATA",),
            rr["C"] + ("Summer16_07Aug2017BCD_V11_DATA",),
            rr["D"] + ("Summer16_07Aug2017BCD_V11_DATA",),
            rr["E"] + ("Summer16_07Aug2017EF_V11_DATA",),
            rr["F"] + ("Summer16_07Aug2017EF_V11_DATA",),
            rr["G"] + ("Summer16_07Aug2017GH_V11_DATA",),
            rr["H"] + ("Summer16_07Aug2017GH_V11_DATA",),
        ],
        "mc": [
            (1, int(1e9), "Summer16_07Aug2017_V11_MC"), # TODO: Or Summer16_23Sep2016V4_MC? (no formula evaluator bug)
        ],
    })

    # https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py
    cfg.set_aux("pileup_mc", [
        1.78653e-05, 2.56602e-05, 5.27857e-05, 8.88954e-05, 0.000109362, 0.000140973,
        0.000240998, 0.00071209, 0.00130121, 0.00245255, 0.00502589, 0.00919534,
        0.0146697, 0.0204126, 0.0267586, 0.0337697, 0.0401478, 0.0450159, 0.0490577,
        0.0524855, 0.0548159, 0.0559937, 0.0554468, 0.0537687, 0.0512055, 0.0476713,
        0.0435312, 0.0393107, 0.0349812, 0.0307413, 0.0272425, 0.0237115, 0.0208329,
        0.0182459, 0.0160712, 0.0142498, 0.012804, 0.011571, 0.010547, 0.00959489,
        0.00891718, 0.00829292, 0.0076195, 0.0069806, 0.0062025, 0.00546581, 0.00484127,
        0.00407168, 0.00337681, 0.00269893, 0.00212473, 0.00160208, 0.00117884,
        0.000859662, 0.000569085, 0.000365431, 0.000243565, 0.00015688, 9.88128e-05,
        6.53783e-05, 3.73924e-05, 2.61382e-05, 2.0307e-05, 1.73032e-05, 1.435e-05,
        1.36486e-05, 1.35555e-05, 1.37491e-05, 1.34255e-05, 1.33987e-05, 1.34061e-05,
        1.34211e-05, 1.34177e-05, 1.32959e-05, 1.33287e-05
    ])

    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II
    cfg.set_aux("min_bias_xs", sn.Number(69.2, (sn.Number.REL, 0.046)))  # mb

    # file merging information (stage -> dataset -> files after merging)
    cfg.set_aux("file_merging", {
        "trees": {
            #"tt_dl": ,
            #"tt_sl": ,
            #"dy_lep_50ToInf": ,
            #"st_tW_t": ,
            #"st_tW_tbar": ,
            #"ttZJets_lep": ,
        }
    })

    # versions
    cfg.set_aux("versions", {
        "WriteTrees": "prod1",
        "MergeTrees": "prod1",
        "MergeMetaData": "prod1",
        "WriteHistograms": "prod1",
        "MergeHistograms": "prod1",
        "MeasureCScaleFactors": "prod1",
        "MeasureScaleFactors": "prod1",
        "FitScaleFactors": "prod1",
        "GetScaleFactorWeights": "prod1",
        "MergeScaleFactorWeights": "prod1",
    })

    return cfg
