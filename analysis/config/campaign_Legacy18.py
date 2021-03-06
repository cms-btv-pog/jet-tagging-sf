# -*- coding: utf-8 -*-
# flake8: noqa

"""
Defintion of the campaign and datasets for 2018 data.
"""


import order as od

from analysis.config.processes import *


# campaign
campaign_name = "Run2_pp_13TeV_Legacy18"
campaign = od.Campaign(
    campaign_name, 3,
    ecm=13,
    bx=25,
)

# datasets
dataset_data_A_ee = od.Dataset(
    "data_A_ee", 0,
    campaign=campaign,
    is_data=True,
    n_files=4790,
    keys=["/EGamma/Run2018A-17Sep2018-v2/MINIAOD"],
    context=campaign_name,
)

dataset_data_B_ee = od.Dataset(
    "data_B_ee", 1,
    campaign=campaign,
    is_data=True,
    n_files=1941,
    keys=["/EGamma/Run2018B-17Sep2018-v1/MINIAOD"],
    context=campaign_name,
)

dataset_data_C_ee = od.Dataset(
    "data_C_ee", 2,
    campaign=campaign,
    is_data=True,
    n_files=2183,
    keys=["/EGamma/Run2018C-17Sep2018-v1/MINIAOD"],
    context=campaign_name,
)

dataset_data_D_ee = od.Dataset(
    "data_D_ee", 3,
    campaign=campaign,
    is_data=True,
    n_files=8510,
    keys=["/EGamma/Run2018D-22Jan2019-v2/MINIAOD"],
    context=campaign_name,
)

datasets_data_ee = [
    dataset_data_A_ee, dataset_data_B_ee, dataset_data_C_ee, dataset_data_D_ee,
]

dataset_data_A_emu = od.Dataset(
    "data_A_emu", 10,
    campaign=campaign,
    is_data=True,
    n_files=541,
    keys=["/MuonEG/Run2018A-17Sep2018-v1/MINIAOD"],
    context=campaign_name,
)

dataset_data_B_emu = od.Dataset(
    "data_B_emu", 11,
    campaign=campaign,
    is_data=True,
    n_files=333,
    keys=["/MuonEG/Run2018B-17Sep2018-v1/MINIAOD"],
    context=campaign_name,
)

dataset_data_C_emu = od.Dataset(
    "data_C_emu", 12,
    campaign=campaign,
    is_data=True,
    n_files=292,
    keys=["/MuonEG/Run2018C-17Sep2018-v1/MINIAOD"],
    context=campaign_name,
)

dataset_data_D_emu = od.Dataset(
    "data_D_emu", 13,
    campaign=campaign,
    is_data=True,
    n_files=1373,
    keys=["/MuonEG/Run2018D-PromptReco-v2/MINIAOD"],
    context=campaign_name,
)

datasets_data_emu = [
    dataset_data_A_emu, dataset_data_B_emu, dataset_data_C_emu, dataset_data_D_emu
]

dataset_data_A_mumu = od.Dataset(
    "data_A_mumu", 20,
    campaign=campaign,
    is_data=True,
    n_files=1260,
    keys=["/DoubleMuon/Run2018A-17Sep2018-v2/MINIAOD"],
    context=campaign_name,
)

dataset_data_B_mumu = od.Dataset(
    "data_B_mumu", 21,
    campaign=campaign,
    is_data=True,
    n_files=614,
    keys=["/DoubleMuon/Run2018B-17Sep2018-v1/MINIAOD"],
    context=campaign_name,
)

dataset_data_C_mumu = od.Dataset(
    "data_C_mumu", 22,
    campaign=campaign,
    is_data=True,
    n_files=526,
    keys=["/DoubleMuon/Run2018C-17Sep2018-v1/MINIAOD"],
    context=campaign_name,
)

dataset_data_D_mumu = od.Dataset(
    "data_D_mumu", 23,
    campaign=campaign,
    is_data=True,
    n_files=2353,
    keys=["/DoubleMuon/Run2018D-PromptReco-v2/MINIAOD"],
    context=campaign_name,
)

datasets_data_mumu = [
    dataset_data_A_mumu, dataset_data_B_mumu, dataset_data_C_mumu, dataset_data_D_mumu,
]

dataset_data_A_mu = od.Dataset(
    "data_A_mu", 30,
    campaign=campaign,
    is_data=True,
    n_files=3224,
    keys=["/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD"],
    context=campaign_name,
)

dataset_data_B_mu = od.Dataset(
    "data_B_mu", 31,
    campaign=campaign,
    is_data=True,
    n_files=1551,
    keys=["/SingleMuon/Run2018B-17Sep2018-v1/MINIAOD"],
    context=campaign_name,
)

dataset_data_C_mu = od.Dataset(
    "data_C_mu", 32,
    campaign=campaign,
    is_data=True,
    n_files=1556,
    keys=["/SingleMuon/Run2018C-17Sep2018-v1/MINIAOD"],
    context=campaign_name,
)

dataset_data_D_mu = od.Dataset(
    "data_D_mu", 33,
    campaign=campaign,
    is_data=True,
    n_files=5195,
    keys=["/SingleMuon/Run2018D-22Jan2019-v2/MINIAOD"],
    context=campaign_name,
)

datasets_data_mu = [
    dataset_data_A_mu, dataset_data_B_mu, dataset_data_C_mu, dataset_data_D_mu,
]

# Single electron data shares the same files with ee
dataset_data_A_e = od.Dataset(
    "data_A_e", 40,
    campaign=campaign,
    is_data=True,
    n_files=dataset_data_A_ee.n_files,
    keys=dataset_data_A_ee.keys,
    context=campaign_name,
)

dataset_data_B_e = od.Dataset(
    "data_B_e", 41,
    campaign=campaign,
    is_data=True,
    n_files=dataset_data_B_ee.n_files,
    keys=dataset_data_B_ee.keys,
    context=campaign_name,
)

dataset_data_C_e = od.Dataset(
    "data_C_e", 42,
    campaign=campaign,
    is_data=True,
    n_files=dataset_data_C_ee.n_files,
    keys=dataset_data_C_ee.keys,
    context=campaign_name,
)

dataset_data_D_e = od.Dataset(
    "data_D_e", 43,
    campaign=campaign,
    is_data=True,
    n_files=dataset_data_D_ee.n_files,
    keys=dataset_data_D_ee.keys,
    context=campaign_name,
)

datasets_data_e = [
    dataset_data_A_e, dataset_data_B_e, dataset_data_C_e, dataset_data_D_e,
]

# MC datasets

# tt

dataset_tt_dl = od.Dataset(
    "tt_dl", 101,
    campaign=campaign,
    n_files=968,
    keys=[
         "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    ],
    context=campaign_name,
)

dataset_tt_sl = od.Dataset(
    "tt_sl", 102,
    campaign=campaign,
    n_files=1523,
    keys=[
        "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    ],
    context=campaign_name,
)

# Drell-Yan

dataset_dy_lep_10To50 = od.Dataset(
    "dy_lep_10To50", 2230,
    campaign=campaign,
    n_files=569,
    keys=[
        "/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    ],
    context=campaign_name,
)

dataset_dy_lep_50ToInf = od.Dataset(
    "dy_lep_50ToInf", 2231,
    campaign=campaign,
    n_files=15 + 2802, # 1254
    keys=[
        "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
        "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM",
        #"/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    ],
    context=campaign_name,
)

dataset_dy_lep_0Jets = od.Dataset(
    "dy_lep_0Jets", 2240,
    campaign=campaign,
    n_files=1405,
    keys=[
        "/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    ],
    context=campaign_name,
)

dataset_dy_lep_1Jets = od.Dataset(
    "dy_lep_1Jets", 2241,
    campaign=campaign,
    n_files=1517,
    keys=[
        "/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    ],
    context=campaign_name,
)

dataset_dy_lep_2Jets = od.Dataset(
    "dy_lep_2Jets", 2242,
    campaign=campaign,
    n_files=1102,
    keys=[
        "/DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    ],
    context=campaign_name,
)

# single top

# s-channel
dataset_st_s_lep = od.Dataset(
    "st_s_lep", 300,
    campaign=campaign,
    n_files=353,
    keys= [
        "/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v4/MINIAODSIM",
    ],
    context=campaign_name,
)

# t-channel
dataset_st_t_t = od.Dataset(
    "st_t_t", 301,
    campaign=campaign,
    n_files=2395,
    keys= [
        "/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    ],
    context=campaign_name,
)

dataset_st_t_tbar = od.Dataset(
    "st_t_tbar", 302,
    campaign=campaign,
    n_files=1319,
    keys= [
        "/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    ],
    context=campaign_name,
)


# tW-channel
dataset_st_tW_t = od.Dataset(
    "st_tW_t", 321,
    campaign=campaign,
    n_files=173,
    keys=[
         "/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM",
    ],
    context=campaign_name,
)

dataset_st_tW_tbar = od.Dataset(
    "st_tW_tbar", 322,
    campaign=campaign,
    n_files=160,
    keys=[
         "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM",
    ],
    context=campaign_name,
)

# diboson

dataset_WW = od.Dataset(
    "WW", 401,
    campaign=campaign,
    n_files=153,
    keys=[
        "/WW_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    ],
    context=campaign_name,
)

dataset_WZ = od.Dataset(
    "WZ", 402,
    campaign=campaign,
    n_files=72,
    keys=[
        "/WZ_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM",
    ],
    context=campaign_name,
)

dataset_ZZ = od.Dataset(
    "ZZ", 403,
    campaign=campaign,
    n_files=53,
    keys=[
        "/ZZ_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    ],
    context=campaign_name,
)

# W + jets

dataset_W_lep = od.Dataset(
    "W_lep", 500,
    campaign=campaign,
    n_files=990,
    keys=[
        "/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    ],
    context=campaign_name,
)

# tt+X

dataset_ttH_bb = od.Dataset(
    "ttH_bb", 601,
    campaign=campaign,
    n_files=319,
    keys=[
        "/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM",
    ],
    context=campaign_name,
)

dataset_ttH_nonbb = od.Dataset(
    "ttH_nonbb", 602,
    campaign=campaign,
    n_files=224,
    keys=[
        "/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    ],
    context=campaign_name,
)

dataset_ttWJets = od.Dataset(
    "ttWJets", 700,
    campaign=campaign,
    n_files=374,
    keys=[
        "/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM",
    ],
    context=campaign_name,
)

dataset_ttZJets = od.Dataset(
    "ttZJets", 710,
    campaign=campaign,
    n_files=632,
    keys=[
        "/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM",
    ],
    context=campaign_name,
)

# link processes to datasets
for d in datasets_data_ee:
    d.add_process(process_data_ee)
for d in datasets_data_emu:
    d.add_process(process_data_emu)
for d in datasets_data_mumu:
    d.add_process(process_data_mumu)
for d in datasets_data_mu:
    d.add_process(process_data_mu)
for d in datasets_data_e:
    d.add_process(process_data_e)

dataset_tt_dl.add_process(process_tt_dl)
dataset_tt_sl.add_process(process_tt_sl)
dataset_dy_lep_10To50.add_process(process_dy_lep_10To50)
dataset_dy_lep_50ToInf.add_process(process_dy_lep_50ToInf)
dataset_dy_lep_0Jets.add_process(process_dy_lep_0Jets)
dataset_dy_lep_1Jets.add_process(process_dy_lep_1Jets)
dataset_dy_lep_2Jets.add_process(process_dy_lep_2Jets)

dataset_st_s_lep.add_process(process_st_s_lep)
dataset_st_t_t.add_process(process_st_t_t)
dataset_st_t_tbar.add_process(process_st_t_tbar)
dataset_st_tW_t.add_process(process_st_tW_t)
dataset_st_tW_tbar.add_process(process_st_tW_tbar)
dataset_WW.add_process(process_WW)
dataset_WZ.add_process(process_WZ)
dataset_ZZ.add_process(process_ZZ)
dataset_W_lep.add_process(process_W_lep)
dataset_ttH_bb.add_process(process_ttH_bb)
dataset_ttH_nonbb.add_process(process_ttH_nonbb)
dataset_ttWJets.add_process(process_ttWJets)
dataset_ttZJets.add_process(process_ttZJets)
