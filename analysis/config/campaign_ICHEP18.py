# -*- coding: utf-8 -*-

"""
Defintion of the campaign and datasets for 2017 data for ICHEP 2018.
"""


import order as od

from analysis.config.processes import process_data_ee, process_data_emu, process_data_mumu, \
    process_tt_dl, process_dy_lep_10To50, process_dy_lep_50ToInf, process_st_tW_t, \
    process_st_tW_tbar, process_WW_sl


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
    n_files=759,
    n_events=58088760,
    keys=["/DoubleEG/Run2017B-17Nov2017-v1/MINIAOD"],
)

dataset_data_C_ee = od.Dataset(
    "data_C_ee", 2,
    campaign=campaign,
    is_data=True,
    n_files=847,
    n_events=65181125,
    keys=["/DoubleEG/Run2017C-17Nov2017-v1/MINIAOD"],
)

dataset_data_D_ee = od.Dataset(
    "data_D_ee", 3,
    campaign=campaign,
    is_data=True,
    n_files=344,
    n_events=25911432,
    keys=["/DoubleEG/Run2017D-17Nov2017-v1/MINIAOD"],
)

dataset_data_E_ee = od.Dataset(
    "data_E_ee", 4,
    campaign=campaign,
    is_data=True,
    n_files=639,
    n_events=56235775,
    keys=["/DoubleEG/Run2017E-17Nov2017-v1/MINIAOD"],
)

dataset_data_F_ee = od.Dataset(
    "data_F_ee", 5,
    campaign=campaign,
    is_data=True,
    n_files=1024,
    n_events=74344288,
    keys=["/DoubleEG/Run2017F-17Nov2017-v1/MINIAOD"],
)

datasets_data_ee = [
    dataset_data_B_ee, dataset_data_C_ee, dataset_data_D_ee, dataset_data_E_ee, dataset_data_F_ee,
]

dataset_data_B_emu = od.Dataset(
    "data_B_emu", 11,
    campaign=campaign,
    is_data=True,
    n_files=70,
    n_events=4453465,
    keys=["/MuonEG/Run2017B-17Nov2017-v1/MINIAOD"],
)

dataset_data_C_emu = od.Dataset(
    "data_C_emu", 12,
    campaign=campaign,
    is_data=True,
    n_files=219,
    n_events=15595214,
    keys=["/MuonEG/Run2017C-17Nov2017-v1/MINIAOD"],
)

dataset_data_D_emu = od.Dataset(
    "data_D_emu", 13,
    campaign=campaign,
    is_data=True,
    n_files=126,
    n_events=9164365,
    keys=["/MuonEG/Run2017D-17Nov2017-v1/MINIAOD"],
)

dataset_data_E_emu = od.Dataset(
    "data_E_emu", 14,
    campaign=campaign,
    is_data=True,
    n_files=239,
    n_events=19043421,
    keys=["/MuonEG/Run2017E-17Nov2017-v1/MINIAOD"],
)

dataset_data_F_emu = od.Dataset(
    "data_F_emu", 15,
    campaign=campaign,
    is_data=True,
    n_files=362,
    n_events=25776363,
    keys=["/MuonEG/Run2017F-17Nov2017-v1/MINIAOD"],
)

datasets_data_emu = [
    dataset_data_B_emu, dataset_data_C_emu, dataset_data_D_emu, dataset_data_E_emu,
    dataset_data_F_emu,
]

dataset_data_B_mumu = od.Dataset(
    "data_B_mumu", 21,
    campaign=campaign,
    is_data=True,
    n_files=164,
    n_events=14501767,
    keys=["/DoubleMuon/Run2017B-17Nov2017-v1/MINIAOD"],
)

dataset_data_C_mumu = od.Dataset(
    "data_C_mumu", 22,
    campaign=campaign,
    is_data=True,
    n_files=466,
    n_events=49636525,
    keys=["/DoubleMuon/Run2017C-17Nov2017-v1/MINIAOD"],
)

dataset_data_D_mumu = od.Dataset(
    "data_D_mumu", 23,
    campaign=campaign,
    is_data=True,
    n_files=242,
    n_events=23075733,
    keys=["/DoubleMuon/Run2017D-17Nov2017-v1/MINIAOD"],
)

dataset_data_E_mumu = od.Dataset(
    "data_E_mumu", 24,
    campaign=campaign,
    is_data=True,
    n_files=596,
    n_events=51589091,
    keys=["/DoubleMuon/Run2017E-17Nov2017-v1/MINIAOD"],
)

dataset_data_F_mumu = od.Dataset(
    "data_F_mumu", 25,
    campaign=campaign,
    is_data=True,
    n_files=875,
    n_events=79756560,
    keys=["/DoubleMuon/Run2017F-17Nov2017-v1/MINIAOD"],
)

datasets_data_mumu = [
    dataset_data_B_mumu, dataset_data_C_mumu, dataset_data_D_mumu, dataset_data_E_mumu,
    dataset_data_F_mumu,
]

dataset_tt_dl = od.Dataset(
    "tt_dl", 101,
    campaign=campaign,
    n_files=164 + 1230,
    n_events=8705576 + 69705626,
    keys=[
        "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/MINIAODSIM",
        "/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM",
    ],
)

dataset_dy_lep_10To50 = od.Dataset(
    "dy_lep_10To50", 221,
    campaign=campaign,
    n_files=508,
    n_events=38832197,
    keys=["/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/MINIAODSIM"],
)

dataset_dy_lep_50ToInf = od.Dataset(
    "dy_lep_50ToInf", 222,
    campaign=campaign,
    n_files=372 + 2752,
    n_events=26923935 + 185998625,
    keys=[
        "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM",
        "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1/MINIAODSIM"
    ],
)

dataset_st_tW_t = od.Dataset(
    "st_tW_t", 321,
    campaign=campaign,
    n_files=167,
    n_events=7636171,
    keys=["/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM"],
)

dataset_st_tW_tbar = od.Dataset(
    "st_tW_tbar", 322,
    campaign=campaign,
    n_files=155,
    n_events=7756300,
    keys=["/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM"],
)

dataset_WW_sl = od.Dataset(
    "WW_sl", 42,
    campaign=campaign,
    n_files=50,
    n_events=1818828,
    keys=["/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM"],
)

# link processes to datasets
for d in datasets_data_ee:
    d.add_process(process_data_ee)
for d in datasets_data_emu:
    d.add_process(process_data_emu)
for d in datasets_data_mumu:
    d.add_process(process_data_mumu)
dataset_tt_dl.add_process(process_tt_dl)
dataset_dy_lep_10To50.add_process(process_dy_lep_10To50)
dataset_dy_lep_50ToInf.add_process(process_dy_lep_50ToInf)
dataset_st_tW_t.add_process(process_st_tW_t)
dataset_st_tW_tbar.add_process(process_st_tW_tbar)
dataset_WW_sl.add_process(process_WW_sl)
