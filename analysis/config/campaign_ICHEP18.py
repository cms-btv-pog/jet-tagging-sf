# -*- coding: utf-8 -*-
# flake8: noqa

"""
Defintion of the campaign and datasets for 2017 data for ICHEP 2018.
"""


import order as od

from analysis.config.processes import *


# campaign
campaign = od.Campaign(
    "2017_Run2_pp_13TeV_ICHEP18", 1,
    ecm=13,
    bx=25,
)

# datasets
dataset_data_B_ee = od.Dataset(
    "data_B_ee", 1,
    campaign=campaign,
    is_data=True,
    n_files=578,
    keys=["/DoubleEG/Run2017B-31Mar2018-v1/MINIAOD"],
)

dataset_data_C_ee = od.Dataset(
    "data_C_ee", 2,
    campaign=campaign,
    is_data=True,
    n_files=682,
    keys=["/DoubleEG/Run2017C-31Mar2018-v1/MINIAOD"],
)

dataset_data_D_ee = od.Dataset(
    "data_D_ee", 3,
    campaign=campaign,
    is_data=True,
    n_files=250,
    keys=["/DoubleEG/Run2017D-31Mar2018-v1/MINIAOD"],
)

dataset_data_E_ee = od.Dataset(
    "data_E_ee", 4,
    campaign=campaign,
    is_data=True,
    n_files=661,
    keys=["/DoubleEG/Run2017E-31Mar2018-v1/MINIAOD"],
)

dataset_data_F_ee = od.Dataset(
    "data_F_ee", 5,
    campaign=campaign,
    is_data=True,
    n_files=867,
    keys=["/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD"],
)

datasets_data_ee = [
    dataset_data_B_ee, dataset_data_C_ee, dataset_data_D_ee, dataset_data_E_ee, dataset_data_F_ee,
]

dataset_data_B_emu = od.Dataset(
    "data_B_emu", 11,
    campaign=campaign,
    is_data=True,
    n_files=48,
    keys=["/MuonEG/Run2017B-31Mar2018-v1/MINIAOD"],
)

dataset_data_C_emu = od.Dataset(
    "data_C_emu", 12,
    campaign=campaign,
    is_data=True,
    n_files=158,
    keys=["/MuonEG/Run2017C-31Mar2018-v1/MINIAOD"],
)

dataset_data_D_emu = od.Dataset(
    "data_D_emu", 13,
    campaign=campaign,
    is_data=True,
    n_files=92,
    keys=["/MuonEG/Run2017D-31Mar2018-v1/MINIAOD"],
)

dataset_data_E_emu = od.Dataset(
    "data_E_emu", 14,
    campaign=campaign,
    is_data=True,
    n_files=242,
    keys=["/MuonEG/Run2017E-31Mar2018-v1/MINIAOD"],
)

dataset_data_F_emu = od.Dataset(
    "data_F_emu", 15,
    campaign=campaign,
    is_data=True,
    n_files=315,
    keys=["/MuonEG/Run2017F-31Mar2018-v1/MINIAOD"],
)

datasets_data_emu = [
    dataset_data_B_emu, dataset_data_C_emu, dataset_data_D_emu, dataset_data_E_emu,
    dataset_data_F_emu,
]

dataset_data_B_mumu = od.Dataset(
    "data_B_mumu", 21,
    campaign=campaign,
    is_data=True,
    n_files=119,
    keys=["/DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD"],
)

dataset_data_C_mumu = od.Dataset(
    "data_C_mumu", 22,
    campaign=campaign,
    is_data=True,
    n_files=463,
    keys=["/DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD"],
)

dataset_data_D_mumu = od.Dataset(
    "data_D_mumu", 23,
    campaign=campaign,
    is_data=True,
    n_files=223,
    keys=["/DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD"],
)

dataset_data_E_mumu = od.Dataset(
    "data_E_mumu", 24,
    campaign=campaign,
    is_data=True,
    n_files=574,
    keys=["/DoubleMuon/Run2017E-31Mar2018-v1/MINIAOD"],
)

dataset_data_F_mumu = od.Dataset(
    "data_F_mumu", 25,
    campaign=campaign,
    is_data=True,
    n_files=847,
    keys=["/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD"],
)

datasets_data_mumu = [
    dataset_data_B_mumu, dataset_data_C_mumu, dataset_data_D_mumu, dataset_data_E_mumu,
    dataset_data_F_mumu,
]

# single electron

dataset_data_B_e = od.Dataset(
    "data_B_e", 31,
    campaign = campaign,
    n_files=499,
    keys=["/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD"],
    is_data=True
)

dataset_data_C_e = od.Dataset(
    "data_C_e", 32,
    campaign = campaign,
    n_files=1179,
    keys=["/SingleElectron/Run2017C-31Mar2018-v1/MINIAOD"],
    is_data=True
)

dataset_data_D_e = od.Dataset(
    "data_D_e", 33,
    campaign = campaign,
    n_files=448,
    keys=["/SingleElectron/Run2017D-31Mar2018-v1/MINIAOD"],
    is_data=True
)

dataset_data_E_e = od.Dataset(
    "data_E_e", 34,
    campaign = campaign,
    n_files=1054,
    keys=["/SingleElectron/Run2017E-31Mar2018-v1/MINIAOD"],
    is_data=True
)

dataset_data_F_e = od.Dataset(
    "data_F_e", 35,
    campaign = campaign,
    n_files=1351,
    keys=["/SingleElectron/Run2017F-31Mar2018-v1/MINIAOD"],
    is_data=True
)

datasets_data_e = [
    dataset_data_B_e, dataset_data_C_e, dataset_data_D_e, dataset_data_E_e,
    dataset_data_F_e,
]

# single muon

dataset_data_B_mu = od.Dataset(
    "data_B_mu", 41,
    campaign = campaign,
    n_files=1059,
    keys=["/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD"],
    is_data=True
)

dataset_data_C_mu = od.Dataset(
    "data_C_mu", 42,
    campaign = campaign,
    n_files=1248,
    keys=["/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD"],
    is_data=True
)

dataset_data_D_mu = od.Dataset(
    "data_D_mu", 43,
    campaign = campaign,
    n_files=607,
    keys=["/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD"],
    is_data=True
)

dataset_data_E_mu = od.Dataset(
    "data_E_mu", 44,
    campaign = campaign,
    n_files=1523,
    keys=["/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD"],
    is_data=True
)

dataset_data_F_mu = od.Dataset(
    "data_F_mu", 45,
    campaign = campaign,
    n_files=2435,
    keys=["/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD"],
    is_data=True
)

datasets_data_mu = [
    dataset_data_B_mu, dataset_data_C_mu, dataset_data_D_mu, dataset_data_E_mu,
    dataset_data_F_mu,
]

# MC datasets

# tt

dataset_tt_dl = od.Dataset(
    "tt_dl", 101,
    campaign=campaign,
    n_files=143 + 1112,
    keys=[
        "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_tt_sl = od.Dataset(
    "tt_sl", 102,
    campaign=campaign,
    n_files=615 + 2176,
    keys=[
        "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
        "/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
    ],
)

# Drell-Yan

# TODO: Add when available (current one is mislabelled)
#dataset_dy_lep_4To50_Ht70To100 = od.Dataset(
#    "dy_lep_4To50_Ht70To100", 2211,
#    campaign=campaign,
#    n_files=121,
#    keys=[
#        "/DYJetsToLL_M-4to50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
#    ],
#)

dataset_dy_lep_4To50_Ht100To200 = od.Dataset(
    "dy_lep_4To50_Ht100To200", 2212,
    campaign=campaign,
    n_files=114 + 15,
    keys=[
        "/DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "/DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_4To50_Ht200To400 = od.Dataset(
    "dy_lep_4To50_Ht200To400", 2213,
    campaign=campaign,
    n_files=19,
    keys=[
        "/DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_4To50_Ht400To600 = od.Dataset(
    "dy_lep_4To50_Ht400To600", 2214,
    campaign=campaign,
    n_files=19,
    keys=[
        "/DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_4To50_Ht600ToInf = od.Dataset(
    "dy_lep_4To50_Ht600ToInf", 2215,
    campaign=campaign,
    n_files=34,
    keys=[
        "/DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_50ToInf_Ht70To100 = od.Dataset(
    "dy_lep_50ToInf_Ht70To100", 2221,
    campaign=campaign,
    n_files=207,
    keys=[
        "/DYJetsToLL_M-50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_50ToInf_Ht100To200 = od.Dataset(
    "dy_lep_50ToInf_Ht100To200", 2222,
    campaign=campaign,
    n_files=153 + 61,
    keys=[
        "/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_50ToInf_Ht200To400 = od.Dataset(
    "dy_lep_50ToInf_Ht200To400", 2223,
    campaign=campaign,
    n_files=202 + 22,
    keys=[
        "/DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "/DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_50ToInf_Ht400To600 = od.Dataset(
    "dy_lep_50ToInf_Ht400To600", 2224,
    campaign=campaign,
    n_files=147 + 23,
    keys=[
        "/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_50ToInf_Ht600To800 = od.Dataset(
    "dy_lep_50ToInf_Ht600To800", 2225,
    campaign=campaign,
    n_files=158,
    keys=[
        "/DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_50ToInf_Ht800To1200 = od.Dataset(
    "dy_lep_50ToInf_Ht800To1200", 2226,
    campaign=campaign,
    n_files=66,
    keys=[
        "/DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_50ToInf_Ht1200To2500 = od.Dataset(
    "dy_lep_50ToInf_Ht1200To2500", 2227,
    campaign=campaign,
    n_files=15,
    keys=[
        "/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_50ToInf_Ht2500ToInf = od.Dataset(
    "dy_lep_50ToInf_Ht2500ToInf", 2228,
    campaign=campaign,
    n_files=10,
    keys=[
        "/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_10To50 = od.Dataset(
    "dy_lep_10To50", 2230,
    campaign=campaign,
    n_files=477,
    keys=[
        "/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_50ToInf = od.Dataset(
    "dy_lep_50ToInf", 2231,
    campaign=campaign,
    n_files=305 + 2126,
    keys=[
        "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
        "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    ],
)

# single top

# s-channel
dataset_st_s_lep = od.Dataset(
    "st_s_lep", 300,
    campaign=campaign,
    n_files=143,
    keys= [
        "/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
    ]
)

# t-channel
dataset_st_t_t = od.Dataset(
    "st_t_t", 301,
    campaign=campaign,
    n_files=111,
    keys= [
        "/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_st_t_tbar = od.Dataset(
    "st_t_tbar", 302,
    campaign=campaign,
    n_files=69,
    keys= [
        "/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
    ],
)


# tW-channel
dataset_st_tW_t = od.Dataset(
    "st_tW_t", 321,
    campaign=campaign,
    n_files=130,
    keys=[
        "/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_st_tW_tbar = od.Dataset(
    "st_tW_tbar", 322,
    campaign=campaign,
    n_files=129,
    keys=[
        "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

# diboson

dataset_WW = od.Dataset(
    "WW", 401,
    campaign=campaign,
    n_files=93,
    keys=[
        "/WW_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_WZ = od.Dataset(
    "WZ", 402,
    campaign=campaign,
    n_files=50,
    keys=[
        "/WZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_ZZ = od.Dataset(
    "ZZ", 403,
    campaign=campaign,
    n_files=23,
    keys=[
        "/ZZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

# W+jets
dataset_W_lep = od.Dataset(
    "W_lep", 500,
    campaign=campaign,
    n_files=514+738,
    keys=[
        "/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
        "/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM",
    ],
)

# tt+X

dataset_ttH = od.Dataset(
    "ttH", 600,
    campaign=campaign,
    n_files=208,
    keys=[
        "/ttH_M125_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_ttWJets_lep = od.Dataset(
    "ttWJets_lep", 700,
    campaign=campaign,
    n_files=90,
    keys=[
        "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_ttWJets_had = od.Dataset(
    "ttWJets_had", 701,
    campaign=campaign,
    n_files=18,
    keys=[
        "/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
    ],
)

dataset_ttZJets_lep = od.Dataset(
    "ttZJets_lep", 710,
    campaign=campaign,
    n_files=130,
    keys=[
        "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    ],
)

dataset_ttZJets_had = od.Dataset(
    "ttZJets_had", 711,
    campaign=campaign,
    n_files=16,
    keys=[
        "/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
    ],
)


# link processes to datasets
for d in datasets_data_ee:
    d.add_process(process_data_ee)
for d in datasets_data_emu:
    d.add_process(process_data_emu)
for d in datasets_data_mumu:
    d.add_process(process_data_mumu)
for d in datasets_data_e:
    d.add_process(process_data_e)
for d in datasets_data_mu:
    d.add_process(process_data_mu)

dataset_tt_dl.add_process(process_tt_dl)
dataset_tt_sl.add_process(process_tt_sl)
#dataset_dy_lep_4To50_Ht70To100.add_process(process_dy_lep_5To50_Ht70To100)  # TODO: 5To50 okay? # TODO: enable
dataset_dy_lep_4To50_Ht100To200.add_process(process_dy_lep_5To50_Ht100To200)  # TODO: 5To50 okay?
dataset_dy_lep_4To50_Ht200To400.add_process(process_dy_lep_5To50_Ht200To400)  # TODO: 5To50 okay?
dataset_dy_lep_4To50_Ht400To600.add_process(process_dy_lep_5To50_Ht400To600)  # TODO: 5To50 okay?
dataset_dy_lep_4To50_Ht600ToInf.add_process(process_dy_lep_5To50_Ht600ToInf)  # TODO: 5To50 okay?
dataset_dy_lep_10To50.add_process(process_dy_lep_10To50)
dataset_dy_lep_50ToInf.add_process(process_dy_lep_50ToInf)
dataset_dy_lep_50ToInf_Ht70To100.add_process(process_dy_lep_50ToInf_Ht70To100)
dataset_dy_lep_50ToInf_Ht100To200.add_process(process_dy_lep_50ToInf_Ht100To200)
dataset_dy_lep_50ToInf_Ht200To400.add_process(process_dy_lep_50ToInf_Ht200To400)
dataset_dy_lep_50ToInf_Ht400To600.add_process(process_dy_lep_50ToInf_Ht400To600)
dataset_dy_lep_50ToInf_Ht600To800.add_process(process_dy_lep_50ToInf_Ht600To800)
dataset_dy_lep_50ToInf_Ht800To1200.add_process(process_dy_lep_50ToInf_Ht800To1200)
dataset_dy_lep_50ToInf_Ht1200To2500.add_process(process_dy_lep_50ToInf_Ht1200To2500)
dataset_dy_lep_50ToInf_Ht2500ToInf.add_process(process_dy_lep_50ToInf_Ht2500ToInf)
dataset_st_s_lep.add_process(process_st_s_lep)
dataset_st_t_t.add_process(process_st_t_t)
dataset_st_t_tbar.add_process(process_st_t_tbar)
dataset_st_tW_t.add_process(process_st_tW_t)
dataset_st_tW_tbar.add_process(process_st_tW_tbar)
dataset_WW.add_process(process_WW)
dataset_WZ.add_process(process_WZ)
dataset_ZZ.add_process(process_ZZ)
dataset_W_lep.add_process(process_W_lep)
dataset_ttH.add_process(process_ttH)
dataset_ttWJets_lep.add_process(process_ttWJets_lep)
dataset_ttWJets_had.add_process(process_ttWJets_had)
dataset_ttZJets_lep.add_process(process_ttZJets_lep)
dataset_ttZJets_had.add_process(process_ttZJets_had)
