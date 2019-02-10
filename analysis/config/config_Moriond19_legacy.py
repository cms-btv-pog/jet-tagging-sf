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
#        ch_emu: [
#            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
#            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
#            "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
#            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
#        ],
        ch_mumu: [
#            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
#            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
#            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
#            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
        ],
        #ch_e: [
        #],
        #ch_mu: [
        #],
    })

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
#
#    # JER
#    cfg.set_aux("jer_version", "Fall17_V3")
#
    # JES
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
#            (1, int(1e9), "Fall17_17Nov2017_V32_MC"),
        ],
    })

#    cfg.set_aux("pileup_mc", [
#        3.39597497605e-05, 6.63688402133e-06, 1.39533611284e-05, 3.64963078209e-05, 6.00872171664e-05,
#        9.33932578027e-05, 0.000120591524486, 0.000128694546198, 0.000361697233219, 0.000361796847553,
#        0.000702474896113, 0.00133766053707, 0.00237817050805, 0.00389825605651, 0.00594546732588,
#        0.00856825906255, 0.0116627396044, 0.0148793350787, 0.0179897368379, 0.0208723871946,
#        0.0232564170641, 0.0249826433945, 0.0262245860346, 0.0272704617569, 0.0283301107549,
#        0.0294006137386, 0.0303026836965, 0.0309692426278, 0.0308818046328, 0.0310566806228,
#        0.0309692426278, 0.0310566806228, 0.0310566806228, 0.0310566806228, 0.0307696426944,
#        0.0300103336052, 0.0288355370103, 0.0273233309106, 0.0264343533951, 0.0255453758796,
#        0.0235877272306, 0.0215627588047, 0.0195825559393, 0.0177296309658, 0.0160560731931,
#        0.0146022004183, 0.0134080690078, 0.0129586991411, 0.0125093292745, 0.0124360740539,
#        0.0123547104433, 0.0123953922486, 0.0124360740539, 0.0124360740539, 0.0123547104433,
#        0.0124360740539, 0.0123387597772, 0.0122414455005, 0.011705203844, 0.0108187105305,
#        0.00963985508986, 0.00827210065136, 0.00683770076341, 0.00545237697118, 0.00420456901556,
#        0.00367513566191, 0.00314570230825, 0.0022917978982, 0.00163221454973, 0.00114065309494,
#        0.000784838366118, 0.000533204105387, 0.000358474034915, 0.000238881117601, 0.0001984254989,
#        0.000157969880198, 0.00010375646169, 6.77366175538e-05, 4.39850477645e-05, 2.84298066026e-05,
#        1.83041729561e-05, 1.17473542058e-05, 7.51982735129e-06, 6.16160108867e-06, 4.80337482605e-06,
#        3.06235473369e-06, 1.94863396999e-06, 1.23726800704e-06, 7.83538083774e-07, 4.94602064224e-07,
#        3.10989480331e-07, 1.94628487765e-07, 1.57888581037e-07, 1.2114867431e-07, 7.49518929908e-08,
#        4.6060444984e-08, 2.81008884326e-08, 1.70121486128e-08, 1.02159894812e-08,
#    ])

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
