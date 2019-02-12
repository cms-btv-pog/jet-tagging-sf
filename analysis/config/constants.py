# -*- coding: utf-8 -*-

"""
Physics constants.
"""


import scinum as sn

# constants
N_LEPS = sn.Number(3)
Z_MASS = sn.Number(91.1876, {"z_mass": 0.0021})
BR_W_HAD = sn.Number(0.6741, {"br_w": 0.0027})
BR_W_LEP = 1 - BR_W_HAD
BR_WW_SL = 2 * BR_W_HAD.mul(BR_W_LEP, rho=-1, inplace=False)
BR_WW_DL = BR_W_LEP**2
BR_WW_FH = BR_W_HAD**2
BR_Z_CLEP = sn.Number(0.033658, {"br_z_clep": 0.000023}) * N_LEPS
BR_H_BB = sn.Number(0.5824, {
    "br_h": ("rel", (0.0065**2+0.0072**2+0.0078**2)**0.5, (0.0065**2+0.0074**2+0.0080**2)**0.5)
})
