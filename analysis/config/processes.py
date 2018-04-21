# -*- coding: utf-8 -*-

"""
Physics processes.
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

process_tt_dl = od.Process(
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
)

process_dy_lep_10To50 = process_dy_lep.add_process(
    "dy_lep_10To50", 221,
    label=r"Drell-Yan, $10 \leq m_{ll} \leq 50$",
    xsecs={
        13: sn.Number(18610),
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

process_WW_sl = process_WW.add_process(
    "WW_sl", 42,
    xsecs={
        13: process_WW.get_xsec(13) * BR_WW_SL,
    },
)
