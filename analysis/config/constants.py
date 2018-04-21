# -*- coding: utf-8 -*-

"""
Physics constants.
"""


import scinum as sn


# constants
N_LEPS = sn.Number(3)
BR_W_HAD = sn.Number(0.6741, {"br_w": 0.0027})
BR_W_LEP = 1 - BR_W_HAD
BR_WW_SL = 2 * BR_W_HAD.mul(BR_W_LEP, rho=-1, inplace=False)
BR_WW_DL = BR_W_LEP**2
BR_WW_FH = BR_W_HAD**2
