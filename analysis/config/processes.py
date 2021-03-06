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

from analysis.config.constants import N_LEPS, BR_W_LEP, BR_Z_CLEP, BR_WW_DL, BR_WW_SL, BR_H_BB


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

process_data_e = od.Process(
    "data_e", 4,
    is_data=True,
    label="data",
)

process_data_mu = od.Process(
    "data_mu", 5,
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

process_tt_sl = process_tt.add_process(
    "tt_sl", 11,
    label=r"$t\bar{t}$ + Jets, SL",
    xsecs={
        13: process_tt.get_xsec(13) * BR_WW_SL,
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

process_dy_lep_50ToInf = process_dy_lep.add_process(
    "dy_lep_50ToInf", 222,
    label=r"Drell-Yan, $m_{ll} > 50$",
    xsecs={
        13: sn.Number(1921.8, {
            "integration": 0.6,
            "pdf": 33.2,
        }) * N_LEPS,
    },
)

process_dy_lep_0Jets = process_dy_lep.add_process(
    "dy_lep_0Jets", 230,
    label=r"Drell-Yan, 0 Jets",
    xsecs={
        13: sn.Number(4620.52),
    },
)

process_dy_lep_1Jets = process_dy_lep.add_process(
    "dy_lep_1Jets", 231,
    label=r"Drell-Yan, 1 Jet",
    xsecs={
        13: sn.Number(859.59),
    },
)

process_dy_lep_2Jets = process_dy_lep.add_process(
    "dy_lep_2Jets", 232,
    label=r"Drell-Yan, 2 Jets",
    xsecs={
        13: sn.Number(338.26),
    },
)

process_st = od.Process(
    "st", 30,
    label=r"Single $t$/$\bar{t}$",
)

# s-channel
process_st_s = process_st.add_process(
    "st_s", 330,
    xsecs={
        13: sn.Number(11.36, {
            "scale": 0.18,
            "pdf": (0.40, 0.45)
        })
    }
)

process_st_s_lep = process_st_s.add_process(
    "st_s_lep", 341,
    xsecs={
        13: process_st_s.get_xsec(13) * BR_W_LEP
    }
)

# t-channel
process_st_t = process_st.add_process(
    "st_t", 310,
    xsecs={
        13: sn.Number(216.99, {
            "scale": (6.62, 4.64),
            "pdf": 6.16,
            "mtop": 1.81
        })
    }
)

process_st_t_t = process_st_t.add_process(
    "st_t_t", 311,
    xsecs={
        13: sn.Number(136.02, {
            "scale": (4.09, 2.92),
            "pdf": 3.52,
            "mtop": 1.11
        })
    }
)

process_st_t_tbar = process_st_t.add_process(
    "st_t_tbar", 312,
    xsecs={
        13: sn.Number(80.95, {
            "scale": (2.53, 1.71),
            "pdf": 3.18,
            "mtop": (0.71, 0.70)
        })
    }
)

# tW-channel
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

# diboson

process_VV = od.Process(
    "VV", 40,
    label=r"VV + Jets",
)

process_WW = process_VV.add_process(
    "WW", 410,
    label="WW",
    xsecs={
        13: sn.Number(118.7, {
            "scale": ("rel", 2.5, 2.2),
            "pdf": ("rel", 0.046),
        }),
    },
)

process_WW_dl = process_WW.add_process(
    "WW_dl", 411,
    xsecs={
        13: process_WW.get_xsec(13) * BR_WW_DL,
    },
)

process_WW_sl = process_WW.add_process(
    "WW_sl", 412,
    xsecs={
        13: process_WW.get_xsec(13) * BR_WW_SL,
    },
)

process_WZ = process_VV.add_process(
    "WZ", 420,
    label="WZ",
    xsecs={
        13: sn.Number.add(
            sn.Number(0.106, {
                "scale": 0.0036,
                "pdf": 0.0050,
            }), sn.Number(0.0663, {
                "scale": 0.0029,
                "pdf": 0.0032,
            }),
            rho=1,
        ) * N_LEPS / BR_W_LEP * N_LEPS / BR_Z_CLEP,
    },
)

process_ZZ = process_VV.add_process(
    "ZZ", 430,
    label="ZZ",
    xsecs={
        13: sn.Number(0.0719, {
            "scale": 0.0023,
            "pdf": 0.0027,
        }) * (N_LEPS / BR_Z_CLEP)**2 / 2.,
    },
)

# W + jets

process_W_lep = od.Process(
    "W_lep", 50,
    label=r"W + Jets, $W \rightarrow l \nu$",
    xsecs={
        13: sn.Number(20508.9, {
            "scale": (165.7, 88.2),
            "pdf": 770.9,
        }) * N_LEPS,
    },
)

# tt + X

process_ttH = od.Process(
    "ttH", 60,
    label=r"$t\bar{t}$ + H",
    xsecs={
        13: sn.Number(0.5071, {
            "scale": ("rel", 0.058, 0.092),
            "pdf": ("rel", 0.036),
        })
    }
)

process_ttH_bb = process_ttH.add_process(
    "ttH_bb", 61,
    xsecs={
        13: process_ttH.get_xsec(13) * BR_H_BB,
    }
)

process_ttH_nonbb = process_ttH.add_process(
    "ttH_nonbb", 62,
    xsecs={
        13: process_ttH.get_xsec(13) * (1 - BR_H_BB),
    }
)

process_ttVJets = od.Process(
    "ttVJets", 70,
    label=r"$t\bar{t}$V + Jets",
)

process_ttWJets = process_ttVJets.add_process(
    "ttWJets", 710,
    label=r"$t\bar{t}$ + W + Jets"
)

process_ttWJets_lep = process_ttWJets.add_process(
    "ttWJets_lep", 711,
    label=r"$t\bar{t}$ + W(lep) + Jets",
    xsecs={
        13: sn.Number(0.2043, {
            "scale": 0.0020
        })
    }
)

process_ttWJets_had = process_ttWJets.add_process(
    "ttWJets_had", 712,
    label=r"$t\bar{t}$ + W(had) + Jets",
    xsecs={
        13: sn.Number(0.4062, {
            "scale": 0.0021
        })
    }
)

process_ttWJets.xsecs = {
    13: sn.Number.add(process_ttWJets_lep.get_xsec(13), process_ttWJets_had.get_xsec(13), rho=1, inplace=False)
}


process_ttZJets = process_ttVJets.add_process(
    "ttZJets", 720,
    label=r"$t\bar{t}$ + Z + Jets"
)

process_ttZJets_lep = process_ttZJets.add_process(
    "ttZJets_lep", 721,
    label=r"$t\bar{t}$ + Z(lep) + Jets",
    xsecs={
        13: sn.Number(0.2529, {
            "scale": 0.0004
        })
    }
)

process_ttZJets_had = process_ttZJets.add_process(
    "ttZJets_had", 722,
    label=r"$t\bar{t}$ + Z(had) + Jets",
    xsecs={
        13: sn.Number(0.5297, {
            "scale": 0.0008
        })
    }
)

process_ttZJets.xsecs = {
    13: sn.Number.add(process_ttZJets_lep.get_xsec(13), process_ttZJets_had.get_xsec(13), rho=1, inplace=False)
}


process_ttVJets.xsecs = {
    13: sn.Number.add(process_ttWJets.get_xsec(13), process_ttZJets.get_xsec(13), rho=1, inplace=False)
}
