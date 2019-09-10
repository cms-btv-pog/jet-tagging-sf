# -*- coding: utf-8 -*-
# flake8: noqa

"""
Definition of the analysis for extracting jet tagging scale factors
"""


import re
import os

import numpy as np
import order as od
import scinum as sn
import six
from order.util import join_root_selection

from analysis.config.constants import Z_MASS
from analysis.config.processes import process_data_e, process_data_mu, process_data_ee, \
    process_data_emu, process_data_mumu, process_tt_dl, process_tt_sl, process_dy_lep, \
    process_st_s, process_st_t, process_st_tW, process_WW, process_WZ, process_ZZ, \
    process_W_lep, process_ttH, process_ttVJets


# define the analysis
analysis = od.Analysis("jet_tagging_sf", 1)

# create base config
cfg = analysis.add_config(name="base", id=0)

# link processes
cfg.add_process(process_data_ee)
cfg.add_process(process_data_emu)
cfg.add_process(process_data_mumu)
cfg.add_process(process_data_e)
cfg.add_process(process_data_mu)
cfg.add_process(process_tt_dl)
cfg.add_process(process_tt_sl)
cfg.add_process(process_dy_lep)
cfg.add_process(process_st_s)
cfg.add_process(process_st_t)
cfg.add_process(process_st_tW)
cfg.add_process(process_WW)
cfg.add_process(process_WZ)
cfg.add_process(process_ZZ)
cfg.add_process(process_W_lep)
cfg.add_process(process_ttH)
cfg.add_process(process_ttVJets)


# define channels
ch_ee = cfg.add_channel("ee", 1)
ch_emu = cfg.add_channel("emu", 2)
ch_mumu = cfg.add_channel("mumu", 3)

# define configurations that are not part of a config

jes_sources = [
    "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation", "SinglePionECAL",
    "SinglePionHCAL", "FlavorQCD", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
    "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal", "RelativeFSR",
    "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF", "PileUpDataMC", "PileUpPtRef",
    "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "Total",
]
jes_total_shifts = {"jesTotal_up", "jesTotal_down"}

cfg.set_aux("jes_sources", jes_sources[:])
# jes sources are required for ShiftTasks already on class level to define list of shift
# however, no information about the config instance is available at that point
if os.environ.get("JTSF_CAMPAIGN", None) is None:
    raise Exception("JTSF campaign has to be defined.")
if os.environ["JTSF_CAMPAIGN"] == "2018_Run2_pp_13TeV_MORIOND19":
    jes_sources.insert(0, "AbsoluteSample")

# add auxiliary info to base config
cfg.set_aux("sandboxes", {
    "slc6": "singularity::/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel6",
    "NO_SANDBOX": "singularity::None",
})

cfg.set_aux("jes_levels", {
    "data": ["L1FastJet", "L2Relative", "L3Absolute", "L2L3Residual"],
    "mc": ["L1FastJet", "L2Relative", "L3Absolute"],
})

cfg.set_aux("btaggers", {
    "deepcsv": {
        "variable": "deepcsv_bcomb",
        "label": "DeepCSV",
    },
    "deepjet": {
        "variable": "deepjet_bcomb",
        "label": "DeepJet",
    },
})

# flavor IDs for .csv result file
cfg.set_aux("flavor_ids", {
    "lf": 2,
    "c": 1,
    "hf": 0,
})

# store binning information
hf_binning = {
    "pt": [20, 30, 50, 70, 100, np.inf],
    "abs(eta)": [0., 2.5],
    "deepcsv": {
        "plotting": [
            -2.01, 0.0, 0.0254, 0.0508, 0.0762, 0.1016, 0.127, 0.1522, 0.2205, 0.2889, 0.3573,
            0.4257, 0.4941, 0.5553, 0.6165, 0.6777, 0.7389, 0.8001, 0.842, 0.884, 0.926, 0.968,
            1.01,
        ],
        "measurement": [
            -2.01, 0.0, 0.1522, 0.2205, 0.2889, 0.3573, 0.4257, 0.4941, 0.5553, 0.6165, 0.6777,
            0.7389, 0.8001, 0.842, 0.884, 0.926, 0.968, 1.01
        ],
    },
    "deepjet": {
        "plotting": [
            -2.01, 0.0, 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.12, 0.16, 0.24,
            0.32, 0.40, 0.48, 0.56, 0.64, 0.72, 0.8, 0.84, 0.88, 0.92, 0.95, 0.98,
            1.01,
        ],
        "measurement": [
            -2.01, 0.0, 0.16, 0.24, 0.32, 0.40, 0.48, 0.56, 0.64, 0.72,
            0.8, 0.84, 0.88, 0.92, 0.95, 0.98, 1.01,
        ],
    },
}

cfg.set_aux("binning", {
    "lf": {
        "pt": [20, 30, 40, 60, np.inf],
        "abs(eta)": [0., 0.8, 1.6, 2.5],
        "deepcsv": {
            "plotting": [
                -2.01, 0.0, 0.0254, 0.0508, 0.0762, 0.1016, 0.127, 0.1522, 0.2205, 0.2889, 0.3573,
                0.4257, 0.4941, 0.5961, 0.6981, 0.8001, 0.835, 0.87, 0.905, 0.94, 0.975, 1.01,
            ],
            "measurement": [
                -2.01, 0.0, 0.0254, 0.0508, 0.0762, 0.1016, 0.127, 0.1522, 0.2205, 0.2889, 0.3573,
                0.4257, 0.4941, 1.01,
            ],
        },
        "deepjet": {
            "plotting": [
                -2.01, 0.0, 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.12, 0.16, 0.24,
                0.32, 0.40, 0.48, 0.56, 0.64, 0.72, 0.8, 0.84, 0.88, 0.92, 0.95, 0.98,
                1.01,
            ],
            "measurement": [
                -2.01, 0.0, 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.12, 0.16, 0.24,
                0.32, 1.01,
            ],
        }
    },
    "hf": hf_binning,
    "c": hf_binning,
})

# information related to systematic shifts
# scaling factors for contamination
cfg.set_aux("contamination_factors", {
    "lf_up": 1.2,
    "lf_down": 0.8,
    "hf_up": 1.2,
    "hf_down": 0.8,
})

# define nested categories (analysis phase space -> hf/lf region -> flavor -> pt bin -> eta bin)
def get_phasespace_info():
    return [
        ("measure", join_root_selection(["n_jets{jec_identifier} == 2", "mll > 12", "dr_ll > 0.2"])),
        ("closure", join_root_selection(["n_jets{jec_identifier} >= 2", "mll > 12", "dr_ll > 0.2"])),
    ]

def get_btag_info(cfg, idx, working_point, b_tagger="deepcsv", operator=">"):
    btag_cfg = cfg.get_aux("btaggers")[b_tagger]
    btag_variable = cfg.get_variable("jet{}_{}".format(idx, btag_cfg["variable"]))
    btag_wp = cfg.get_aux("working_points")[b_tagger][working_point]

    return "({}) {} {}".format(btag_variable.expression, operator, btag_wp)

def get_z_window_info(cfg, flavour, et_miss=30.0, z_window=10.):
    cuts = []
    # ET-miss requirement
    et_miss_expr = "(met_px{jec_identifier}**2 + met_py{jec_identifier}**2)**0.5"
    if flavour == "hf":
        cuts.append("{} > {}".format(et_miss_expr, et_miss))
    elif flavour == "lf":
        cuts.append("{} < {}".format(et_miss_expr, et_miss))
    else:
        raise ValueError("Unrecognized flavour {}".format(flavour))

    # z-mass window
    if flavour == "hf":
        cuts.append("abs(mll - {}) > {}".format(Z_MASS.nominal, z_window))
    else:
        cuts.append("abs(mll - {}) < {}".format(Z_MASS.nominal, z_window))

    # z peak diamond
    if flavour == "lf":
        cuts.append("pass_z_mask{jec_identifier} == 0")
    return cuts

def get_region_info(cfg, idx, channel, et_miss=30., z_window=10., add_btag_cut=True, b_tagger="deepcsv"):
    hf_cuts, lf_cuts = [], []
    if add_btag_cut:
        hf_cuts.append(get_btag_info(cfg, idx, "medium", b_tagger, ">"))
        lf_cuts.append(get_btag_info(cfg, idx, "loose", b_tagger, "<"))

    if channel != "emu":
        lf_cuts.extend(get_z_window_info(cfg, "lf", et_miss=et_miss, z_window=z_window))
        hf_cuts.extend(get_z_window_info(cfg, "hf", et_miss=et_miss, z_window=z_window))

    return [
        ("hf", join_root_selection(hf_cuts)),
        ("lf", join_root_selection(lf_cuts)),
    ]

def get_contamination_region_info(cfg, channel, et_miss=30.0, z_window=10.0, b_tagger="deepcsv"):
    cuts = []
    cuts.append(get_btag_info(cfg, 1, "tight", b_tagger, ">"))
    cuts.append(get_btag_info(cfg, 2, "tight", b_tagger, ">"))

    if channel != "emu":
        cuts.extend(get_z_window_info(cfg, "lf", et_miss=et_miss, z_window=z_window))

    return [
        ("cont", join_root_selection(cuts)),
    ]

def get_flavor_info(idx):
    return [
        ("b", "abs(jet%d_flavor{jec_identifier}) == 5" % idx),
        ("c", "abs(jet%d_flavor{jec_identifier}) == 4" % idx),
        ("udsg", "abs(jet%d_flavor{jec_identifier}) != 5 && abs(jet%d_flavor{jec_identifier}) != 4" % (idx, idx)),
        ("inclusive", "1 * 1"),
    ]

def binning_to_selection(binning, variable):
    def to_string(value):
        if np.isinf(value):
            return "Inf"
        elif value == 0:
            return "0"
        else:
            return str(value).replace(".", "p")

    selections = []
    for left_edge, right_edge in zip(binning[:-1], binning[1:]):
        name = "{}To{}".format(to_string(left_edge), to_string(right_edge))
        cuts = []
        if left_edge != 0:
            cuts.append("{} > {}".format(variable, left_edge))
        if not np.isinf(right_edge):
            cuts.append("{} <= {}".format(variable, right_edge))
        edges = [left_edge, right_edge]
        selections.append((name, join_root_selection(cuts), edges))
    return selections

def get_axis_info(cfg, idx, axis_var, fmt=None):
    if fmt is None:
        fmt = "jet{}_" + axis_var

    binning = cfg.get_aux("binning")
    hf_bins = binning["hf"][axis_var]
    lf_bins = binning["lf"][axis_var]
    variable = fmt.format(idx)
    return {
        "hf": binning_to_selection(hf_bins, variable),
        "lf": binning_to_selection(lf_bins, variable),
    }

def get_category(cfg, pt, eta, region, b_tagger, phase_space="measure"):
    matches = []
    for category in cfg.categories:
        if not category.has_tag(b_tagger):
            continue
        cat_phasespace = category.get_aux("phase_space", None)
        if not cat_phasespace == phase_space:
            continue

        cat_region = category.get_aux("region", None)
        if not region == cat_region:
            continue

        cat_pt_range = category.get_aux("pt", (0., 0.))
        if not (cat_pt_range[0] < pt <= cat_pt_range[1]):
            continue

        cat_eta_range = category.get_aux("eta", (0., 0.))
        if not (cat_eta_range[0] < eta <= cat_eta_range[1]):
            continue
        matches.append(category)
    # only return the category if the matching is unambiguos
    if len(matches) == 1:
        return matches[0]
    else:
        raise ValueError("Expected one single matching category, but got {}".format(matches))

# variables
cfg.add_variable(
    name="dr_ll",
    expression="dr_ll",
    binning=(25, 0., 5.,),
    x_title="dR(ll)",
)
cfg.add_variable(
    name="mll",
    expression="mll",
    binning=(20, 80., 100.),
    tags={"contamination"},
    x_title="M(ll)",
)
cfg.add_variable(
    name="n_jets",
    expression="n_jets",
    binning=(10, 0., 10.,),
    x_title="N(jets)",
)
cfg.add_variable(
    name="n_tags_deepcsv",
    expression="n_tags_deepcsv",
    aux={"b_tagger": "deepcsv"}, # to filter required b-tagger in histogram writer
    binning=(10, 0., 10.,),
    tags={"n_tags"},
)
cfg.add_variable(
    name="n_tags_deepjet",
    expression="n_tags_deepjet",
    aux={"b_tagger": "deepjet"}, # to filter required b-tagger in histogram writer
    binning=(10, 0., 10.,),
    tags={"n_tags"},
)

for lep_idx in xrange(1, 3):
    cfg.add_variable(
        name="lep{}_pt".format(lep_idx),
        expression="(lep{}_px**2 + lep{}_py**2)**0.5".format(lep_idx, lep_idx),
        binning=(25, 0., 500.,),
        unit="GeV",
        tags={"main"},
        x_title="Lep_{{{}}} p_{{T}}".format(lep_idx),
    )

for jet_idx in xrange(1, 5):
    tags = {"basic"}
    if jet_idx <= 2:
        tags = tags | {"main"}
    cfg.add_variable(
        name="jet{}_pt".format(jet_idx),
        expression="jet{}_pt{{jec_identifier}}".format(jet_idx),
        binning=(25, 0., 500.,),
        unit="GeV",
        x_title="Jet_{{{}}} p_{{T}}".format(jet_idx),
        tags=tags
    )
    cfg.add_variable(
        name="jet{}_eta".format(jet_idx),
        expression="jet{}_eta{{jec_identifier}}".format(jet_idx),
        binning=(25, -2.5, 2.5),
        x_title="Jet_{{{}}} Eta".format(jet_idx),
        tags=tags,
    )


def add_btag_variables(cfg):
    for jet_idx in xrange(1, 5):
        for region in [None, "hf", "lf"]:
            if not region:
                binning = (25, 0., 1.)
                tags = {"skip_lf", "skip_hf", "basic"}
                postfix = ""
            elif region == "hf":
                binning = cfg.get_aux("binning")["hf"]["deepcsv"]["plotting"]
                tags = {"skip_lf", "basic"}
                postfix = "_hf"
            elif region == "lf":
                binning = cfg.get_aux("binning")["lf"]["deepcsv"]["plotting"]
                tags = {"skip_hf", "basic"}
                postfix = "_lf"

            tags = tags | {"b_tagging"}
            if jet_idx <= 2:
                tags = tags | {"measurement", "main"}
            cfg.add_variable(
                name="jet{}_deepcsv_b{}".format(jet_idx, postfix),
                expression="jet{}_deepcsv_b{{jec_identifier}}".format(jet_idx),
                binning=binning,
                x_title="Jet_{{{}}} prob_{{b}}".format(jet_idx),
                context=cfg.name,
            )
            cfg.add_variable(
                name="jet{}_deepcsv_bb{}".format(jet_idx, postfix),
                expression="jet{}_deepcsv_bb{{jec_identifier}}".format(jet_idx),
                binning=binning,
                x_title="Jet_{{{}}} prob_{{bb}}".format(jet_idx),
                context=cfg.name,
            )
            # deepcsv discriminator
            cfg.add_variable(
                name="jet{}_deepcsv_bcomb{}".format(jet_idx, postfix),
                expression="jet{0}_deepcsv_b{{jec_identifier}} + jet{0}_deepcsv_bb{{jec_identifier}}".format(jet_idx),
                binning=binning,
                x_title="Jet_{{{}}} DeepCSV".format(jet_idx),
                tags=tags,
                aux={"b_tagger": "deepcsv"}, # to filter required b-tagger in histogram writer
                context=cfg.name,
            )

            # deepjet discriminator
            if region in ["lf", "hf"]:
                binning = cfg.get_aux("binning")[region]["deepjet"]["plotting"]
            cfg.add_variable(
                name="jet{}_deepjet_bcomb{}".format(jet_idx, postfix),
                expression="jet{0}_deepjet_b{{jec_identifier}} + jet{0}_deepjet_bb{{jec_identifier}} + "\
                    "jet{0}_deepjet_lepb{{jec_identifier}}".format(jet_idx),
                binning=binning,
                x_title="Jet_{{{}}} DeepJet".format(jet_idx),
                tags=tags,
                aux={"b_tagger": "deepjet"}, # to filter required b-tagger in histogram writer
                context=cfg.name,
            )

def add_categories(cfg, b_tagger):
    # categories
    for ch in [ch_ee, ch_emu, ch_mumu]:
        # phase space region loop (measurement, closure, ...)
        for ps_name, ps_sel in get_phasespace_info():
            # inclusive region categories to measure rates
            for rg_name, rg_sel in get_region_info(cfg, 1, ch, add_btag_cut=False, b_tagger=b_tagger):
                # we skip the emu channel in the lf region because the DY (the main contribution)
                # should have same-flavored leptons
                if rg_name == "lf" and ch == ch_emu:
                    continue

                # categories to perform overall normalization of each channel
                rg_cat_combined = ch.add_category(
                    name="{}__{}__{}__{}__{}".format(ch.name, ps_name, rg_name, b_tagger, cfg.name),
                    label="{}, {}, {}".format(ch.name, ps_name, rg_name),
                    selection=join_root_selection("channel == {}".format(ch.id), ps_sel, rg_sel),
                    tags={"scales", b_tagger},
                    aux={
                        "channel": ch,
                        "phase_space": ps_name,
                        "region": rg_name,
                        "config": cfg.name,
                    },
                )
                # combine region categories to create inclusive control regions for plotting
                rg_merged_name = "{}__{}__{}".format(ps_name, rg_name, b_tagger)
                if not cfg.has_category(rg_merged_name):
                    rg_merged_cat = cfg.add_category(
                        name=rg_merged_name,
                        label="{}, {}".format(ps_name, rg_name),
                        tags={"inclusive", b_tagger},
                        aux={
                            "phase_space": ps_name,
                            "region": rg_name,
                        },
                        context=cfg.name,
                    )
                else:
                    rg_merged_cat = cfg.get_category(rg_merged_name)
                rg_merged_cat.add_category(rg_cat_combined)

            # loop over both jet1 jet2 permutations
            for i_tag_jet, i_probe_jet in [(1, 2), (2, 1)]:
                # region loop (hf, lf, ...)
                for rg_name, rg_sel in get_region_info(cfg, i_tag_jet, ch, b_tagger=b_tagger):
                    if rg_name == "lf" and ch == ch_emu:
                        continue

                    rg_cat = ch.add_category(
                        name="{}__{}__{}__j{}__{}__{}".format(ch.name, ps_name, rg_name, i_tag_jet, b_tagger, cfg.name),
                        label="{}, {}, {} region (j{} tagged)".format(ch.name, ps_name, rg_name, i_tag_jet),
                        selection=join_root_selection("channel == {}".format(ch.id), ps_sel, rg_sel),
                        tags={b_tagger},
                    )

                    # combined region categories, with tag jet cut applied
                    # used to determine e.g. sample composition in measurement regions
                    rg_btag_merged_name = "{}__{}__{}__{}__btag".format(ps_name, rg_name, b_tagger, cfg.name)
                    if not cfg.has_category(rg_btag_merged_name):
                        rg_btag_merged_cat = cfg.add_category(
                            name=rg_btag_merged_name,
                            label="{}, {}".format(ps_name, rg_name, b_tagger),
                            tags={"combined", b_tagger},
                            aux={
                                "phase_space": ps_name,
                                "region": rg_name,
                            },
                            context=cfg.name,
                        )
                    else:
                        rg_btag_merged_cat = cfg.get_category(rg_btag_merged_name)

                    # flavor loop (b, c, udsg, ...)
                    for fl_name, fl_sel in get_flavor_info(i_probe_jet):
                        fl_cat = rg_cat.add_category(
                            name="{}__f{}".format(rg_cat.name, fl_name),
                            label="{}, {} flavor".format(rg_cat.label, fl_name),
                            selection=join_root_selection(rg_cat.selection, fl_sel),
                            tags={b_tagger},
                        )

                        # pt loop
                        for pt_idx, (pt_name, pt_sel, pt_range) in enumerate(
                            get_axis_info(cfg, i_probe_jet, "pt", "jet{}_pt{{jec_identifier}}")[rg_name]):

                            pt_cat = fl_cat.add_category(
                                name="{}__pt{}".format(fl_cat.name, pt_name),
                                label="{}, pt {}".format(fl_cat.label, pt_name),
                                selection=join_root_selection(fl_cat.selection, pt_sel),
                                tags={b_tagger},
                            )

                            # eta loop
                            for eta_idx, (eta_name, eta_sel, eta_range) in enumerate(
                                get_axis_info(cfg, i_probe_jet, "abs(eta)", fmt="abs(jet{}_eta{{jec_identifier}})")[rg_name]):

                                eta_cat = pt_cat.add_category(
                                    name="{}__eta{}".format(pt_cat.name, eta_name),
                                    label="{}, eta {}".format(pt_cat.label, eta_name),
                                    selection=join_root_selection(pt_cat.selection, eta_sel),
                                    aux={
                                        "channel": ch,
                                        "i_probe_jet": i_probe_jet,
                                        "i_tag_jet": i_tag_jet,
                                        "phase_space": ps_name,
                                        "region": rg_name,
                                        "flavor": fl_name,
                                        "config": cfg.name,
                                    },
                                    tags={b_tagger},
                                )

                                # merged category for both jets and all flavors
                                merged_vars = (ps_name, rg_name, pt_name, eta_name, b_tagger)
                                merged_name = "{}__{}__pt{}__eta{}__{}".format(*merged_vars)

                                # define categories for testing
                                merged_tags = {"merged", b_tagger}
                                if rg_name == "hf" and (pt_idx == 1 and eta_idx == 0):
                                    merged_tags = merged_tags | {"test"}
                                if rg_name == "lf" and (pt_idx == 2 and eta_idx == 0):
                                    merged_tags = merged_tags | {"test"}

                                if not cfg.has_category(merged_name):
                                    label = "{}, {} region, pt {}, eta {}".format(*merged_vars)
                                    merged_cat = cfg.add_category(
                                        name=merged_name,
                                        label=label,
                                        tags=merged_tags,
                                        aux={
                                            "phase_space": ps_name,
                                            "region": rg_name,
                                            "eta": eta_range,
                                            "pt": pt_range,
                                        },
                                        context=cfg.name,
                                    )
                                    if rg_name == "hf":
                                        # add c categories (not written to histograms)
                                        c_vars = (ps_name, "c", pt_name, eta_name, b_tagger)
                                        c_name = "{}__{}__pt{}__eta{}__{}".format(*c_vars)
                                        label = "{}, {} region, pt {}, eta {}".format(*c_vars)
                                        c_cat = cfg.add_category(
                                            name=c_name,
                                            label=label,
                                            tags={"c", b_tagger},
                                            aux={
                                                "phase_space": ps_name,
                                                "region": "c",
                                                "eta": eta_range,
                                                "pt": pt_range,
                                            },
                                            context=cfg.name,
                                        )
                                        c_cat.set_aux("binning_category", merged_cat)

                                else:
                                    merged_cat = cfg.get_category(merged_name)
                                merged_cat.add_category(eta_cat)
                                rg_btag_merged_cat.add_category(eta_cat)

                                # Specialized b-tag discriminant binnings are defined on
                                # the merged categories, but needed when writing leaf categories
                                eta_cat.set_aux("binning_category", merged_cat)

            # add categories to measure light flavour contamination uncertainty
            for rg_name, rg_sel in get_contamination_region_info(cfg, ch, b_tagger=b_tagger):
                if ch == ch_emu:
                    continue

                contamination_cat = ch.add_category(
                    name="{}__{}__{}__{}__{}".format(ch.name, ps_name, rg_name, b_tagger, cfg.name),
                    label="{}, {}, {}".format(ch.name, ps_name, rg_name),
                    selection=join_root_selection("channel == {}".format(ch.id), ps_sel, rg_sel),
                    tags={b_tagger},
                    aux={
                        "channel": ch,
                        "phase_space": ps_name,
                        "region": rg_name,
                        "config": cfg.name,
                    },
                )
                # combine contamination regions over all channels
                cont_merged_name = "{}__{}__{}".format(ps_name, rg_name, b_tagger)
                if not cfg.has_category(cont_merged_name):
                    cont_merged_cat = cfg.add_category(
                        name=cont_merged_name,
                        label="{}, {}".format(ps_name, rg_name),
                        tags={"contamination", b_tagger},
                        aux={
                            "phase_space": ps_name,
                            "region": rg_name,
                        },
                        context=cfg.name,
                    )
                else:
                    cont_merged_cat = cfg.get_category(cont_merged_name)
                cont_merged_cat.add_category(contamination_cat)

def get_file_merging(cfg, key, dataset):
    dataset_name = dataset if isinstance(dataset, six.string_types) else dataset.name
    return cfg.get_aux("file_merging")[key].get(dataset_name, 1)

cfg.set_aux("get_file_merging", get_file_merging)

# add specific configs
from analysis.config.config_ICHEP18 import create_config as create_config_ICHEP18
config_ICHEP18 = create_config_ICHEP18(cfg)
add_btag_variables(config_ICHEP18)
add_categories(config_ICHEP18, "deepcsv")
add_categories(config_ICHEP18, "deepjet")

from analysis.config.config_Moriond19 import create_config as create_config_Moriond19
config_Moriond19 = create_config_Moriond19(cfg)
add_btag_variables(config_Moriond19)
add_categories(config_Moriond19, "deepcsv")
add_categories(config_Moriond19, "deepjet")
config_Moriond19.get_aux("jes_sources").insert(0, "AbsoluteSample")

from analysis.config.config_Moriond19_legacy import create_config as create_config_Moriond19_legacy
config_Moriond19_legacy = create_config_Moriond19_legacy(cfg)
add_btag_variables(config_Moriond19_legacy)
add_categories(config_Moriond19_legacy, "deepcsv")
add_categories(config_Moriond19_legacy, "deepjet")
config_Moriond19_legacy.get_aux("binning")["lf"]["abs(eta)"] = [0., 0.8, 1.6, 2.4]
config_Moriond19_legacy.get_aux("binning")["hf"]["abs(eta)"] = [0., 2.4]
config_Moriond19_legacy.get_aux("binning")["c"]["abs(eta)"] = [0., 2.4]
