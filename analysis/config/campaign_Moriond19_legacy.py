# -*- coding: utf-8 -*-
# flake8: noqa

"""
Defintion of the campaign and datasets for 2017 data for ICHEP 2018.
"""


import order as od

from analysis.config.processes import *


# campaign
campaign = od.Campaign(
    "2018_Run2_pp_13TeV_MORIOND19legacy", 1,
    ecm=13,
    bx=25,
)

# datasets

dataset_data_B_ee = od.Dataset(
    "data_B_ee", 1,
    campaign=campaign,
    is_data=True,
    n_files=922,
    keys=["/DoubleEG/Run2016B-17Jul2018_ver2-v1/MINIAOD"],
)

dataset_data_C_ee = od.Dataset(
    "data_C_ee", 2,
    campaign=campaign,
    is_data=True,
    n_files=427,
    keys=["/DoubleEG/Run2016C-17Jul2018-v1/MINIAOD"],
)

dataset_data_D_ee = od.Dataset(
    "data_D_ee", 3,
    campaign=campaign,
    is_data=True,
    n_files=471,
    keys=["/DoubleEG/Run2016D-17Jul2018-v1/MINIAOD"],
)

dataset_data_E_ee = od.Dataset(
    "data_E_ee", 4,
    campaign=campaign,
    is_data=True,
    n_files=375,
    keys=["/DoubleEG/Run2016E-17Jul2018-v1/MINIAOD"],
)

dataset_data_F_ee = od.Dataset(
    "data_F_ee", 5,
    campaign=campaign,
    is_data=True,
    n_files=309,
    keys=["/DoubleEG/Run2016F-17Jul2018-v1/MINIAOD"],
)

dataset_data_G_ee = od.Dataset(
    "data_G_ee", 6,
    campaign=campaign,
    is_data=True,
    n_files=715,
    keys=["/DoubleEG/Run2016G-17Jul2018-v1/MINIAOD"],
)

dataset_data_H_ee = od.Dataset(
    "data_H_ee", 7,
    campaign=campaign,
    is_data=True,
    n_files=736,
    keys=["/DoubleEG/Run2016H-17Jul2018-v1/MINIAOD"],
)

datasets_data_ee = [
    dataset_data_B_ee, dataset_data_C_ee, dataset_data_D_ee, dataset_data_E_ee,
    dataset_data_F_ee, dataset_data_G_ee, dataset_data_H_ee
]

dataset_data_B_emu = od.Dataset(
    "data_B_emu", 11,
    campaign=campaign,
    is_data=True,
    n_files=249,
    keys=["/MuonEG/Run2016B-17Jul2018_ver2-v1/MINIAOD"],
)

dataset_data_C_emu = od.Dataset(
    "data_C_emu", 12,
    campaign=campaign,
    is_data=True,
    n_files=112,
    keys=["/MuonEG/Run2016C-17Jul2018-v1/MINIAOD"],
)

dataset_data_D_emu = od.Dataset(
    "data_D_emu", 13,
    campaign=campaign,
    is_data=True,
    n_files=192,
    keys=["/MuonEG/Run2016D-17Jul2018-v1/MINIAOD"],
)

dataset_data_E_emu = od.Dataset(
    "data_E_emu", 14,
    campaign=campaign,
    is_data=True,
    n_files=209,
    keys=["/MuonEG/Run2016E-17Jul2018-v2/MINIAOD"],
)

dataset_data_F_emu = od.Dataset(
    "data_F_emu", 15,
    campaign=campaign,
    is_data=True,
    n_files=159,
    keys=["/MuonEG/Run2016F-17Jul2018-v1/MINIAOD"],
)

dataset_data_G_emu = od.Dataset(
    "data_G_emu", 16,
    campaign=campaign,
    is_data=True,
    n_files=302,
    keys=["/MuonEG/Run2016G-17Jul2018-v1/MINIAOD"],
)

dataset_data_H_emu = od.Dataset(
    "data_H_emu", 17,
    campaign=campaign,
    is_data=True,
    n_files=267,
    keys=["/MuonEG/Run2016H-17Jul2018-v1/MINIAOD"],
)

datasets_data_emu = [
    dataset_data_B_emu, dataset_data_C_emu, dataset_data_D_emu, dataset_data_E_emu,
    dataset_data_F_emu, dataset_data_G_emu, dataset_data_H_emu
]

dataset_data_B_mumu = od.Dataset(
    "data_B_mumu", 21,
    campaign=campaign,
    is_data=True,
    n_files=451,
    keys=["/DoubleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD"],
)

dataset_data_C_mumu = od.Dataset(
    "data_C_mumu", 22,
    campaign=campaign,
    is_data=True,
    n_files=203,
    keys=["/DoubleMuon/Run2016C-17Jul2018-v1/MINIAOD"],
)

dataset_data_D_mumu = od.Dataset(
    "data_D_mumu", 23,
    campaign=campaign,
    is_data=True,
    n_files=215,
    keys=["/DoubleMuon/Run2016D-17Jul2018-v1/MINIAOD"],
)

dataset_data_E_mumu = od.Dataset(
    "data_E_mumu", 24,
    campaign=campaign,
    is_data=True,
    n_files=186,
    keys=["/DoubleMuon/Run2016E-17Jul2018-v1/MINIAOD"],
)

dataset_data_F_mumu = od.Dataset(
    "data_F_mumu", 25,
    campaign=campaign,
    is_data=True,
    n_files=155,
    keys=["/DoubleMuon/Run2016F-17Jul2018-v1/MINIAOD"],
)

dataset_data_G_mumu = od.Dataset(
    "data_G_mumu", 26,
    campaign=campaign,
    is_data=True,
    n_files=346,
    keys=["/DoubleMuon/Run2016G-17Jul2018-v1/MINIAOD"],
)

dataset_data_H_mumu = od.Dataset(
    "data_H_mumu", 27,
    campaign=campaign,
    is_data=True,
    n_files=378,
    keys=["/DoubleMuon/Run2016H-17Jul2018-v1/MINIAOD"],
)

datasets_data_mumu = [
    dataset_data_B_mumu, dataset_data_C_mumu, dataset_data_D_mumu, dataset_data_E_mumu,
    dataset_data_F_mumu, dataset_data_G_mumu, dataset_data_H_mumu
]

# MC datasets

# tt

dataset_tt_dl = od.Dataset(
    "tt_dl", 101,
    campaign=campaign,
    n_files=592,
    keys=[
        "/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
    ],
)

dataset_tt_sl = od.Dataset(
    "tt_sl", 102,
    campaign=campaign,
    n_files=1142,
    keys=[
        "/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
    ],
)

# Drell-Yan

dataset_dy_lep_10To50 = od.Dataset(
    "dy_lep_10To50", 2230,
    campaign=campaign,
    n_files=264,
    keys=[
        "/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
    ],
)

dataset_dy_lep_50ToInf = od.Dataset(
    "dy_lep_50ToInf", 2231,
    campaign=campaign,
    n_files=360+701,
    keys=[
        "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM",
        "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM",
    ],
)

# single top

# s-channel
dataset_st_s_lep = od.Dataset(
    "st_s_lep", 300,
    campaign=campaign,
    n_files=8,
    keys=[
        "/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
    ]
)

# t-channel
dataset_st_t_t = od.Dataset(
    "st_t_t", 301,
    campaign=campaign,
    n_files=499,
    keys= [
        "/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
    ],
)

dataset_st_t_tbar = od.Dataset(
    "st_t_tbar", 302,
    campaign=campaign,
    n_files=298,
    keys= [
        "/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM",
    ],
)


# tW-channel
dataset_st_tW_t = od.Dataset(
    "st_tW_t", 321,
    campaign=campaign,
    n_files=60,
    keys=[
        "/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM",
    ],
)

dataset_st_tW_tbar = od.Dataset(
    "st_tW_tbar", 322,
    campaign=campaign,
    n_files=58,
    keys=[
        "/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM",
    ],
)

# diboson

dataset_WW = od.Dataset(
    "WW", 401,
    campaign=campaign,
    n_files=7+53,
    keys=[
        "/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM",
    ],
)

dataset_WZ = od.Dataset(
    "WZ", 402,
    campaign=campaign,
    n_files=8+29,
    keys=[
        "/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM",
    ],
)

dataset_ZZ = od.Dataset(
    "ZZ", 403,
    campaign=campaign,
    n_files=7,
    keys=[
        "/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
    ],
)

# W + jets

dataset_W_lep = od.Dataset(
    "W_lep", 500,
    campaign=campaign,
    n_files=215+410,
    keys=[
        "/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
        "/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/MINIAODSIM"
    ],
)

# tt+X

# TODO: ttH bb+nonbb

dataset_ttWJets_lep = od.Dataset(
    "ttWJets_lep", 700,
    campaign=campaign,
    n_files=31,
    keys=[
        "/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM",
    ],
)

dataset_ttWJets_had = od.Dataset(
    "ttWJets_had", 701,
    campaign=campaign,
    n_files=7,
    keys=[
        "/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
    ],
)

dataset_ttZJets_lep = od.Dataset(
    "ttZJets_lep", 710,
    campaign=campaign,
    n_files=49+48,
    keys=[
        "/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM",
        "/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext3-v1/MINIAODSIM",
    ],
)

dataset_ttZJets_had = od.Dataset(
    "ttZJets_had", 711,
    campaign=campaign,
    n_files=7,
    keys=[
        "/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
    ],
)

# link processes to datasets
for d in datasets_data_ee:
    d.add_process(process_data_ee)
for d in datasets_data_emu:
    d.add_process(process_data_emu)
for d in datasets_data_mumu:
    d.add_process(process_data_mumu)


dataset_tt_dl.add_process(process_tt_dl)
dataset_tt_sl.add_process(process_tt_sl)
dataset_dy_lep_10To50.add_process(process_dy_lep_10To50)
dataset_dy_lep_50ToInf.add_process(process_dy_lep_50ToInf)
dataset_st_s_lep.add_process(process_st_s_lep)
dataset_st_t_t.add_process(process_st_t_t)
dataset_st_t_tbar.add_process(process_st_t_tbar)
dataset_st_tW_t.add_process(process_st_tW_t)
dataset_st_tW_tbar.add_process(process_st_tW_tbar)
dataset_WW.add_process(process_WW)
dataset_WZ.add_process(process_WZ)
dataset_ZZ.add_process(process_ZZ)
dataset_W_lep.add_process(process_W_lep)
dataset_ttWJets_lep.add_process(process_ttWJets_lep)
dataset_ttWJets_had.add_process(process_ttWJets_had)
dataset_ttZJets_lep.add_process(process_ttZJets_lep)
dataset_ttZJets_had.add_process(process_ttZJets_had)
