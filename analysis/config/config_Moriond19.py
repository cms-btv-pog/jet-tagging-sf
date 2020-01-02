# -*- coding: utf-8 -*-

import scinum as sn

def create_config(base_cfg):
    # setup the config for Moriond 2019 (2018 data)
    from analysis.config.campaign_Moriond19 import campaign as campaign_Moriond19
    from analysis.config.jet_tagging_sf import ch_ee, ch_emu, ch_mumu
    cfg = base_cfg.copy(campaign=campaign_Moriond19)
    ch_e = cfg.add_channel("e", 4)
    ch_mu = cfg.add_channel("mu", 5)

    # add datasets
    dataset_names = [
        "data_A_ee", "data_B_ee", "data_C_ee", "data_D_ee",
        "data_A_emu", "data_B_emu", "data_C_emu", "data_D_emu",
        "data_A_mumu", "data_B_mumu", "data_C_mumu", "data_D_mumu",
        "data_A_e", "data_B_e", "data_C_e", "data_D_e",
        "data_A_mu", "data_B_mu", "data_C_mu", "data_D_mu",
        "tt_dl", "tt_sl",
        #"dy_lep_10To50",
        #"dy_lep_50ToInf",
        "dy_lep_0Jets", "dy_lep_1Jets", "dy_lep_2Jets",
        "st_s_lep",
        "st_t_t", "st_t_tbar",
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
        ch_ee: 59966.1613198,
        ch_emu: 59966.1613198,
        ch_mumu: 59966.1613198,
        ch_e: 59966.1613198,
        ch_mu: 59966.1613198,
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
        "data": "102X_dataRun2_v12", # 102X_dataRun2_Prompt_v15
        "mc": "102X_upgrade2018_realistic_v20",
    })

    # lumi, normtag and pileup file
    cfg.set_aux("lumi_file", "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/"
        "ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt")
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
    # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2?rev=131
    cfg.set_aux("metFilters", {
        "data": [
            "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter",
            "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter", #"Flag_BadChargedCandidateFilter",
            "Flag_eeBadScFilter", "Flag_ecalBadCalibReducedMINIAODFilter",
        ],
        "mc": [
            "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter",
            "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter", #"Flag_BadChargedCandidateFilter",
            "Flag_ecalBadCalibReducedMINIAODFilter",
        ],
    })

    # JER
    # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    cfg.set_aux("jer_version", "Autumn18_V7")

    # JES
    # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
    cfg.set_aux("jes_version", {
        "data": [
            rr["A"] + ("Autumn18_RunA_V19_DATA",),
            rr["B"] + ("Autumn18_RunB_V19_DATA",),
            rr["C"] + ("Autumn18_RunC_V19_DATA",),
            rr["D"] + ("Autumn18_RunD_V19_DATA",),
        ],
        "mc": [
            (1, int(1e9), "Autumn18_V19_MC"),
        ],
    })

    cfg.set_aux("jes_uncertainty_file", {
        "factorized": None,  # take file from jes github
        "reduced": "https://cernbox.cern.ch/index.php/s/abTn1z0lpGZugkn/download",
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
            "data_D_mumu": 2,
            "data_D_e": 3,
            "data_A_mu": 2,
            "data_D_mu": 3,
            "tt_dl": 484,
            "tt_sl": 508,
            "dy_lep_50ToInf": 149,
            "dy_lep_0Jets": 14,
            "dy_lep_1Jets": 85,
            "dy_lep_2Jets": 221,
            "st_s_lep": 21,
            "st_t_t": 55,
            "st_t_tbar": 30,
            "st_tW_t": 22,
            "st_tW_tbar": 18,
            "WW": 3,
            "WZ": 2,
            "ZZ": 2,
            "W_lep": 5,
            "ttH_bb": 46,
            "ttH_nonbb": 38,
            "ttWJets": 63,
            "ttZJets": 106
        }
    })

    # versions
    cfg.set_aux("versions", {
        "WriteTrees": "prod15", # reduced jes sources
        "MergeTrees": "prod15",
        "MergeMetaData": "prod15",
        "WriteHistograms": "prod15",
        "MergeHistograms": "prod15",
        "MeasureCScaleFactors": "prod15",
        "MeasureScaleFactors": "prod15",
        "FitScaleFactors": "prod15",
        "BundleScaleFactors": "prod15",
        "GetScaleFactorWeights": "prod15",
        "MergeScaleFactorWeights": "prod15",
        "OptimizeBinning": "prod1",
        "CreateScaleFactorResults": "prod15",
    })

    # add sl categories
    def add_categories(cfg, b_tagger):
        from order.util import join_root_selection
        from analysis.config.jet_tagging_sf import get_flavor_info

        # categories
        for ch in [ch_e, ch_mu]:
            # phase space region loop (measurement, closure, ...)
            for ps_name, ps_sel in [("closure", join_root_selection(["n_jets{jec_identifier} == 4",
                    "n_tags_{}{{jec_identifier}} == 2".format(b_tagger)]))]:
                for jet_idx in range(1, 5):
                    for fl_name, fl_sel in get_flavor_info(jet_idx):
                        # categories per channel
                        rg_cat = ch.add_category(
                            name="{}__{}__j{}__{}__{}__{}".format(ch.name, ps_name, str(jet_idx),
                                fl_name, b_tagger, cfg.name),
                            label="{}, {}, jet{}, {}".format(ch.name, ps_name, str(jet_idx), fl_name),
                            selection=join_root_selection("channel == {}".format(ch.id), ps_sel, fl_sel),
                            tags={b_tagger},
                            aux={
                                "channel": ch,
                                "phase_space": ps_name,
                                "config": cfg.name,
                                "flavor": fl_name,
                                "i_flavor_jet": jet_idx,
                            },
                        )
                        # combine region categories to create inclusive control regions for plotting
                        rg_merged_name = "sl__{}__{}".format(ps_name, b_tagger)
                        if not cfg.has_category(rg_merged_name):
                            rg_merged_cat = cfg.add_category(
                                name=rg_merged_name,
                                label="sl, {}".format(ps_name),
                                tags={"sl", b_tagger},
                                aux={
                                    "phase_space": ps_name,
                                },
                                context=cfg.name,
                            )
                        else:
                            rg_merged_cat = cfg.get_category(rg_merged_name)
                        rg_merged_cat.add_category(rg_cat)

    add_categories(cfg, "deepcsv")
    add_categories(cfg, "deepjet")
    return cfg
