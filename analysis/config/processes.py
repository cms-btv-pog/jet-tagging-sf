# -*- coding: utf-8 -*-

"""
Physics processes.
If not stated otherwise, cross sections are given in pb.
Values taken from:
- https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat13TeVInclusive?rev=18
- https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns?rev=151
"""


import order as od
import scinum as sn

from analysis.config.constants import N_LEPS, BR_WW_DL, BR_WW_SL


process_data_ee = od.Process(
    "data_ee", 1,
    is_data=True,
    label="data",
)

process_data_emu = od.Process(
    "data_emu", 2,
    is_data=True,
    label="data",
)

process_data_mumu = od.Process(
    "data_mumu", 3,
    is_data=True,
    label="data",
)

process_tt = od.Process(
    "tt", 10,
    label=r"$t\bar{t}$ + Jets",
    xsecs={
        13: sn.Number(831.76, {
            "scale": (19.77, 29.20),
            "pdf": 35.06,
            "mtop": (23.18, 22.45),
        }),
    },
)

process_tt_dl = process_tt.add_process(
    "tt_dl", 12,
    label=r"$t\bar{t}$ + Jets, DL",
    xsecs={
        13: process_tt.get_xsec(13) * BR_WW_DL,
    },
)

process_dy = od.Process(
    "dy", 20,
    label="Drell-Yan",
)

process_dy_lep = process_dy.add_process(
    "dy_lep", 22,
    label=r"Drell-Yan, $Z \rightarrow ll$",
)

process_dy_lep_5To50 = process_dy_lep.add_process(
    "dy_lep_5To50", 221,
    label=r"Drell-Yan, $5 \leq m_{ll} \leq 50$",
    xsecs={
        13: sn.Number(71310, {
            "scale": 70,
        }),
    },

)

process_dy_lep_10To50 = process_dy_lep.add_process(
    "dy_lep_10To50", 223,
    label=r"Drell-Yan, $10 \leq m_{ll} \leq 50$",
    xsecs={
        13: sn.Number(18610),
    },
)

process_dy_lep_5To50_Ht70To100 = process_dy_lep_5To50.add_process(
    "dy_lep_5To50_Ht70To100", 2211,
    xsecs={
        13: sn.Number(301.2, {
            "pdf": 0.8,
        }),
    },
)

process_dy_lep_5To50_Ht100To200 = process_dy_lep_5To50.add_process(
    "dy_lep_5To50_Ht100To200", 2212,
    xsecs={
        13: sn.Number(224.2, {
            "pdf": 5.7,
        }),
    },
)

process_dy_lep_5To50_Ht200To400 = process_dy_lep_5To50.add_process(
    "dy_lep_5To50_Ht200To400", 2213,
    xsecs={
        13: sn.Number(37.2, {
            "pdf": 1.1,
        }),
    },
)

process_dy_lep_5To50_Ht400To600 = process_dy_lep_5To50.add_process(
    "dy_lep_5To50_Ht400To600", 2214,
    xsecs={
        13: sn.Number(3.581, {
            "pdf": 0.118,
        }),
    },
)

process_dy_lep_5To50_Ht600ToInf = process_dy_lep_5To50.add_process(
    "dy_lep_5To50_Ht600ToInf", 2215,
    xsecs={
        13: sn.Number(1.124, {
            "pdf": 0.038,
        }),
    },
)

process_dy_lep_50ToInf = process_dy_lep.add_process(
    "dy_lep_50ToInf", 222,
    label=r"Drell-Yan, $m_{ll} \gt 50$",
    xsecs={
        13: sn.Number(1921.8, {
            "integration": 0.6,
            "pdf": 33.2,
        }) * N_LEPS,
    },
)

process_dy_lep_50ToInf_Ht70To100 = process_dy_lep_50ToInf.add_process(
    "dy_lep_50ToInf_Ht70To100", 2221,
    xsecs={
        13: sn.Number(169.9, {
            "pdf": 0.5,
        }),
    },
)

process_dy_lep_50ToInf_Ht100To200 = process_dy_lep_50ToInf.add_process(
    "dy_lep_50ToInf_Ht100To200", 2222,
    xsecs={
        13: sn.Number(147.40, {
            "pdf": 0.09,
        }),
    },
)

process_dy_lep_50ToInf_Ht200To400 = process_dy_lep_50ToInf.add_process(
    "dy_lep_50ToInf_Ht200To400", 2223,
    xsecs={
        13: sn.Number(40.99, {
            "pdf": 0.04,
        }),
    },
)

process_dy_lep_50ToInf_Ht400To600 = process_dy_lep_50ToInf.add_process(
    "dy_lep_50ToInf_Ht400To600", 2224,
    xsecs={
        13: sn.Number(5.678, {
            "pdf": 0.005,
        }),
    },
)

process_dy_lep_50ToInf_Ht600To800 = process_dy_lep_50ToInf.add_process(
    "dy_lep_50ToInf_Ht600To800", 2225,
    xsecs={
        13: sn.Number(1.367),
    },
)

process_dy_lep_50ToInf_Ht800To1200 = process_dy_lep_50ToInf.add_process(
    "dy_lep_50ToInf_Ht800To1200", 2226,
    xsecs={
        13: sn.Number(0.6304),
    },
)

process_dy_lep_50ToInf_Ht1200To2500 = process_dy_lep_50ToInf.add_process(
    "dy_lep_50ToInf_Ht1200To2500", 2227,
    xsecs={
        13: sn.Number(0.1514),
    },
)

process_dy_lep_50ToInf_Ht2500ToInf = process_dy_lep_50ToInf.add_process(
    "dy_lep_50ToInf_Ht2500ToInf", 2228,
    xsecs={
        13: sn.Number(0.003565),
    },
)

process_st = od.Process(
    "st", 30,
    label=r"Single $t$/$\bar{t}$",
)

process_st_tW = process_st.add_process(
    "st_tW", 320,
    label=r"Single $t$/$\bar{t}$, tW-channel",
    xsecs={
        13: sn.Number(71.7, {
            "scale": 1.8,
            "pdf": 3.4,
        }),
    },
)

process_st_tW_t = process_st_tW.add_process(
    "st_tW_t", 321,
    xsecs={
        13: sn.Number(35.85, {
            "scale": 0.90,
            "pdf": 1.70,
        }),
    },
)

process_st_tW_tbar = process_st_tW.add_process(
    "st_tW_tbar", 322,
    xsecs={
        13: sn.Number(35.85, {
            "scale": 0.90,
            "pdf": 1.70,
        }),
    },
)

process_WW = od.Process(
    "WW", 40,
    label="WW",
    xsecs={
        13: sn.Number(118.7, {
            "scale": ("rel", 2.5, 2.2),
            "pdf": ("rel", 0.046),
        }),
    },
)

process_WW_dl = process_WW.add_process(
    "WW_dl", 41,
    xsecs={
        13: process_WW.get_xsec(13) * BR_WW_DL,
    },
)

process_WW_sl = process_WW.add_process(
    "WW_sl", 42,
    xsecs={
        13: process_WW.get_xsec(13) * BR_WW_SL,
    },
)
