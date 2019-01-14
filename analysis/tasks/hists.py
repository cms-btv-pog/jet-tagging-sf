# -*- coding: utf-8 -*-


import os
import collections
import itertools
import array

import law
import luigi
import six
import numpy as np
from order.util import join_root_selection
from collections import defaultdict

from analysis.config.jet_tagging_sf import get_category, jes_sources
from analysis.tasks.base import AnalysisTask, DatasetTask, ShiftTask, WrapperTask, GridWorkflow
from analysis.tasks.trees import MergeTrees, MergeMetaData
from analysis.tasks.external import CalculatePileupWeights
from analysis.util import TreeExtender


class WriteHistograms(DatasetTask, GridWorkflow, law.LocalWorkflow):

    iteration = luigi.IntParameter(default=0, description="iteration of the scale factor "
        "calculation, starting at zero, default: 0")
    final_it = luigi.BoolParameter(description="Flag for the final iteration of the scale factor "
        "calculation.")
    variable_tag = luigi.Parameter(default=None, description="Only consider variables with the given "
        "tag. Use all if empty.")
    shifts = CSVParameter(default=[])

    file_merging = "trees"

    workflow_run_decorators = [law.decorator.notify]

    def __init__(self, *args, **kwargs):
        super(WriteHistograms, self).__init__(*args, **kwargs)
        # set shifts
        if self.dataset_inst.is_data:
            shifts = {"nominal"}
        else:
            shifts = {"nominal"} | {"jes{}_{}".format(shift, direction) for shift, direction in itertools.product(
                jes_sources, ["up", "down"])}
            if self.iteration > 0:
                shifts = shifts | {"{}_{}".format(shift, direction) for shift, direction in itertools.product(
                    ["lf", "hf", "lf_stats1", "lf_stats2", "hf_stats1", "hf_stats2"], ["up", "down"])}
        if len(self.shifts) == 0:
            self.shifts = shifts
        elif any([shift not in shifts for shift in self.shifts]):
            raise ValueError("Unknown shift in {}".format(self.shifts))

    def workflow_requires(self):
        from analysis.tasks.measurement import FitScaleFactors

        reqs = super(WriteHistograms, self).workflow_requires()

        if not self.cancel_jobs and not self.cleanup_jobs:
            reqs["meta"] = MergeMetaData.req(self, version=self.get_version(MergeMetaData),
                _prefer_cli=["version"])
            if self.dataset_inst.is_mc:
                reqs["pu"] = CalculatePileupWeights.req(self)
            if not self.pilot:
                reqs["tree"] = MergeTrees.req(self, cascade_tree=-1,
                    version=self.get_version(MergeTrees), _prefer_cli=["version"])
            if self.iteration > 0:
                reqs["sf"] = {shift: FitScaleFactors.req(self, iteration=self.iteration - 1,
                    shift=shift, fix_normalization=self.final_it,
                    version=self.get_version(FitScaleFactors), _prefer_cli=["version"])
                    for shift in self.shifts}

        return reqs

    def requires(self):
        from analysis.tasks.measurement import FitScaleFactors

        reqs = {
            "tree": MergeTrees.req(self, cascade_tree=self.branch, branch=0,
                version=self.get_version(MergeTrees), _prefer_cli=["version", "workflow"]),
            "meta": MergeMetaData.req(self, version=self.get_version(MergeMetaData),
                _prefer_cli=["version"]),
        }
        if self.dataset_inst.is_mc:
            reqs["pu"] = CalculatePileupWeights.req(self)
        if self.iteration > 0:
            reqs["sf"] = {shift: FitScaleFactors.req(self, iteration=self.iteration - 1,
                shift=shift, fix_normalization=self.final_it,
                version=self.get_version(FitScaleFactors), _prefer_cli=["version"])
                for shift in self.shifts}
        return reqs

    def store_parts(self):
        return super(WriteHistograms, self).store_parts() + (self.iteration,)

    def output(self):
        return self.wlcg_target("hists_{}.root".format(self.branch))

    def get_jec_identifier(self, shift):
        if shift.startswith("jes"):
            return "_" + shift
        else:
            return ""

    def get_pileup_weighter(self, inp, meta):
        with inp.load() as pu_file: # TODO: Generalize, use either mc pileup or calculate
            pu_hist = pu_file.Get("pileup_data")
            pu_values_data = [
                pu_hist.GetBinContent(i)
                for i in range(1, pu_hist.GetNbinsX() + 1)
            ]
        pu_values_mc = meta.load()["pileup"]
        pu_values_mc = [count / sum(pu_values_mc) for count in pu_values_mc]
        pu_values = [ pu_frac_data / pu_frac_mc if pu_frac_mc > 0 else 0. for pu_frac_mc, pu_frac_data in zip(pu_values_mc, pu_values_data)]

        def add_branch(extender):
            extender.add_branch("pu_weight", unpack="pu")

        def add_value(entry):
            # some events have inf pileup, skip them
            weight = 1.
            pu = entry.pu[0]
            if np.isfinite(pu):
                pu_idx = int(pu) - 1
                if 0 <= pu_idx < len(pu_values):
                    weight = pu_values[pu_idx]
            entry.pu_weight[0] = weight

        return add_branch, add_value

    def get_scale_factor_weighter(self, inp, shift):
        with inp.load() as sfs:
            sf_hists = {}
            for category in sfs.GetListOfKeys():
                category_dir = sfs.Get(category.GetName())
                hist = category_dir.Get("sf")
                # decouple from open file
                hist.SetDirectory(0)

                sf_hists[category.GetName()] = hist

        btag_var = self.config_inst.get_aux("btagger")["variable"]
        identifier = self.get_jec_identifier(shift)

        def add_branch(extender):
            unpack_vars = sum(
                [["jet{}_pt{}".format(idx, identifier), "jet{}_flavor{}".format(idx, identifier),
                "jet{}_eta{}".format(idx, identifier), "jet{}_{}{}".format(idx, btag_var, identifier)]
                for idx in range(1, 5)],
                []
            )
            extender.add_branch("scale_factor_lf_{}".format(shift), unpack=unpack_vars)
            extender.add_branch("scale_factor_c_{}".format(shift), unpack=unpack_vars)
            extender.add_branch("scale_factor_hf_{}".format(shift), unpack=unpack_vars)

        def add_value(entry):
            scale_factor_lf = 1.
            scale_factor_c = 1.
            scale_factor_hf = 1.
            for jet_idx in range(1, 5):
                jet_pt = getattr(entry, "jet{}_pt{}".format(jet_idx, identifier))[0]
                jet_eta = getattr(entry, "jet{}_eta{}".format(jet_idx, identifier))[0]
                jet_flavor = getattr(entry, "jet{}_flavor{}".format(jet_idx, identifier))[0]
                jet_btag = getattr(entry, "jet{}_{}{}".format(jet_idx, btag_var, identifier))[0]

                # stop when number of jets is exceeded
                if jet_flavor < -999.:
                    break

                # find category in which the scale factor of the jet was computed to get correct histogram
                # TODO: Handle c-jets
                region = "hf" if abs(jet_flavor) in (4, 5) else "lf"
                category = get_category(jet_pt, abs(jet_eta), region, phase_space="measure")

                # get scale factor
                sf_hist = sf_hists[category.name]
                bin_idx = sf_hist.FindBin(jet_btag)
                scale_factor = sf_hist.GetBinContent(bin_idx)

                if abs(jet_flavor) == 5:
                    scale_factor_hf *= scale_factor
                elif abs(jet_flavor) == 4:
                    scale_factor_c *= 1.  # TODO: set to scale_factor once we have sf hists for c jets
                else:
                    scale_factor_lf *= scale_factor

            getattr(entry, "scale_factor_lf_{}".format(shift))[0] = scale_factor_lf
            getattr(entry, "scale_factor_c_{}".format(shift))[0] = scale_factor_c
            getattr(entry, "scale_factor_hf_{}".format(shift))[0] = scale_factor_hf

        return add_branch, add_value

    @law.decorator.notify
    def run(self):
        import ROOT

        inp = self.input()
        outp = self.output()
        outp.parent.touch(0o0770)

        # get child categories
        categories = []
        channels = [self.config_inst.get_aux("dataset_channels")[self.dataset_inst]] \
            if self.dataset_inst.is_data else self.config_inst.channels.values()
        for channel in channels:
            for category, _, children in channel.walk_categories():
                if not children:
                    categories.append((channel, category))
        categories = list(set(categories))

        # get processes
        if len(self.dataset_inst.processes) != 1:
            raise NotImplementedError("only datasets with exactly one linked process can be"
                " handled, got {}".format(len(self.dataset_inst.processes)))
        processes = list(self.dataset_inst.processes.values())

        # build a progress callback
        progress = self.create_progress_callback(len(categories))

        # open the output file
        with outp.localize("w") as tmp:
            with tmp.dump("RECREATE") as output_file:
                with self.publish_step("creating root output file directories ..."):
                    process_dirs = {}
                    for _, category in categories:
                        output_file.cd()
                        category_dir = output_file.mkdir(category.name)
                        for process in processes:
                            category_dir.cd()
                            process_dir = category_dir.mkdir(process.name)
                            process_dir.Write()
                            process_dirs[(category.name, process.name)] = process_dir

                # open the input file and get the tree
                # as we need to extend the tree with custom weights, we do not cache the file
                with inp["tree"].load("UPDATE", cache=False) as input_file:
                    tree = input_file.Get("tree")
                    self.publish_message("{} events in tree".format(tree.GetEntries()))

                    # identifier for jec shifted variables
                    for shift in self.shifts:
                        jec_identifier = self.get_jec_identifier(shift)

                        # pt aliases for jets
                        for obj in ["jet1", "jet2", "jet3", "jet4"]:
                            tree.SetAlias("{0}_pt{1}".format(obj, jec_identifier),
                                "({0}_px{1}**2 + {0}_py{1}**2)**0.5".format(obj, jec_identifier))
                        # b-tagging alias
                        btag_var = self.config_inst.get_aux("btagger")["variable"]
                        for obj in ["jet1", "jet2", "jet3", "jet4"]:
                            variable = self.config_inst.get_variable("{0}_{1}".format(obj, btag_var))
                            tree.SetAlias(variable.name + jec_identifier, variable.expression.format(
                                **{"jec_identifier": jec_identifier}))
                    # pt aliases for leptons
                    for obj in ["lep1", "lep2"]:
                        tree.SetAlias("{0}_pt".format(obj),
                            "({0}_px**2 + {0}_py**2)**0.5".format(obj))

                    # extend the tree
                    if self.dataset_inst.is_mc:
                        with self.publish_step("extending the input tree with weights ..."):
                            weighters = []

                            # pileup weight
                            weighters.append(self.get_pileup_weighter(inp["pu"], inp["meta"]))

                            # weights from previous iterations
                            if self.iteration > 0:
                                # b-tagging scale factors
                                for shift in self.shifts:
                                    weighters.append(self.get_scale_factor_weighter(
                                        inp["sf"][shift]["sf"], shift))

                            input_file.cd()
                            with TreeExtender(tree) as te:
                                for add_branch, _ in weighters:
                                    add_branch(te)
                                for entry in te:
                                    for _, add_value in weighters:
                                        add_value(entry)

                    for i, (channel, category) in enumerate(categories):
                        self.publish_message("writing histograms in category {} ({}/{})".format(
                            category.name, i + 1, len(categories)))

                        # get the region (HF / LF)
                        # not all child categories have regions associated, e.g. the phase space
                        # inclusive regions ("measure", "closure")
                        region = category.get_aux("region", None)

                        # set weights that are common for all shifts
                        base_weights = []
                        if self.dataset_inst.is_mc: # TODO: Should be done separately for each process
                            base_weights.append("gen_weight")
                            # lumi weight
                            lumi = self.config_inst.get_aux("lumi")[channel]
                            x_sec = process.get_xsec(self.config_inst.campaign.ecm).nominal
                            sum_weights = inp["meta"].load()["event_weights"]["sum"]
                            lumi_weight = lumi * x_sec / sum_weights
                            base_weights.append(str(lumi_weight))

                            # pu weight
                            base_weights.append("pu_weight")

                        for process in processes:
                            # change into the correct directory
                            process_dirs[(category.name, process.name)].cd()
                            for shift in self.shifts:
                                jec_identifier = self.get_jec_identifier(shift)

                                # weights
                                weights = base_weights[:]
                                if self.dataset_inst.is_mc:
                                    # channel scale weight
                                    if self.iteration > 0:
                                        # b-tag scale factor weights
                                        phase_space = category.get_aux("phase_space", None)
                                        # In measurement categories,
                                        # apply scale factors only for contamination
                                        if phase_space == "measure" and not self.final_it:
                                            weights.append("scale_factor_c_{}".format(shift))
                                            if region == "hf":
                                                weights.append("scale_factor_lf_{}".format(shift))
                                            elif region == "lf":
                                                weights.append("scale_factor_hf_{}".format(shift))
                                            else:
                                                raise ValueError("Unexpected region {}".format(region))
                                        else:
                                            weights.append("scale_factor_lf_{}".format(shift))
                                            weights.append("scale_factor_c_{}".format(shift))
                                            weights.append("scale_factor_hf_{}".format(shift))

                                # totalWeight alias
                                while len(weights) < 2:
                                    weights.insert(0, "1")
                                tree.SetAlias("totalWeight", join_root_selection(weights, op="*"))

                                # actual projecting
                                for variable in self.config_inst.variables:
                                    if variable.has_tag("skip_all"):
                                        continue
                                    if region and variable.has_tag("skip_{}".format(region)):
                                        continue
                                    # if a vaiable tag is given, require it
                                    if self.variable_tag is not None and not variable.has_tag(self.variable_tag):
                                        continue

                                    hist = ROOT.TH1F("{}_{}".format(variable.name, shift),
                                        variable.full_title(root=True), variable.n_bins,
                                        array.array("f", variable.bin_edges))
                                    hist.Sumw2()

                                    # build the full selection string, including the total event weight
                                    selection = [
                                        category.selection,
                                        "jetmet_pass{jec_identifier} == 1",
                                        "{} != -10000".format(variable.expression),
                                    ]
                                    if variable.selection:
                                        selection.append(variable.selection)
                                    selection = join_root_selection(selection).format(
                                        **{"jec_identifier": jec_identifier})
                                    selection = join_root_selection(selection, "totalWeight", op="*")

                                    # project and write the histogram
                                    tree.Project("{}_{}".format(variable.name, shift),
                                        variable.expression.format(**{"jec_identifier": jec_identifier}),
                                        selection)
                                    hist.Write()

                        progress(i)


class WriteHistogramsWrapper(WrapperTask):

    wrapped_task = WriteHistograms


class MergeHistograms(AnalysisTask, law.CascadeMerge):

    iteration = WriteHistograms.iteration
    final_it = WriteHistograms.final_it

    merge_factor = 12

    def create_branch_map(self):
        return law.CascadeMerge.create_branch_map(self)

    def cascade_workflow_requires(self):
        reqs = collections.OrderedDict()

        for dataset in self.config_inst.datasets:
            reqs[dataset.name] = WriteHistograms.req(self, dataset=dataset.name,
                version=self.get_version(WriteHistograms), _prefer_cli=["version", "workflow"])

        return reqs

    def trace_cascade_workflow_inputs(self, inputs):
        targets = []
        for d in inputs.values():
            targets.extend(list(d["collection"].targets.values()))
        return targets

    def cascade_requires(self, start_leaf, end_leaf):
        slices = []
        for dataset in self.config_inst.datasets:
            file_merging = WriteHistograms.file_merging
            n_files = self.config_inst.get_aux("get_file_merging")(file_merging, dataset)
            slices.extend([(dataset, i) for i in range(n_files)])

        target_slices = slices[start_leaf:end_leaf]
        return [
            WriteHistograms.req(self, dataset=d.name, branch=l,
                version=self.get_version(WriteHistograms), _prefer_cli=["version"])
            for d, l in target_slices
        ]

    def store_parts(self):
        return super(MergeHistograms, self).store_parts() + (self.iteration,)

    def cascade_output(self):
        return self.wlcg_target("hists.root")

    def merge(self, inputs, output):
        tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
        tmp_dir.touch()

        with output.localize("w") as tmp_out:
            # fetch inputs
            with self.publish_step("fetching inputs ..."):
                def fetch(inp):
                    inp.copy_to_local(tmp_dir.child(inp.unique_basename, type="f"), cache=False)
                    return inp.unique_basename

                def callback(i):
                    self.publish_message("fetch file {} / {}".format(i + 1, len(inputs)))

                bases = law.util.map_verbose(fetch, inputs, every=5, callback=callback)

            with self.publish_step("merging ..."):
                if len(bases) == 1:
                    self.publish_message("only 1 file to merge")
                    tmp_out.path = tmp_dir.child(bases[0]).path
                else:
                    # merge using hadd
                    bases = " ".join(bases)
                    cmd = "hadd -n 0 -d {} {} {}".format(tmp_dir.path, tmp_out.path, bases)
                    code = law.util.interruptable_popen(cmd, shell="True", executable="/bin/bash",
                        cwd=tmp_dir.path)[0]
                    if code != 0:
                        raise Exception("hadd failed")

                    self.publish_message("merged file size: {:.2f} {}".format(
                        *law.util.human_bytes(os.stat(tmp_out.path).st_size)))

    def glite_output_postfix(self):
        return "_{}_{}".format(self.cascade_tree, self.cascade_depth)

    def arc_output_postfix(self):
        return self.glite_output_postfix()


# TODO: Remove redundancies between GetScaleFactorWeights and WriteHistograms?
class GetScaleFactorWeights(DatasetTask, GridWorkflow, law.LocalWorkflow):

    iteration = WriteHistograms.iteration
    file_merging = WriteHistograms.file_merging

    normalize_cerrs = luigi.BoolParameter()

    def __init__(self, *args, **kwargs):
        super(GetScaleFactorWeights, self).__init__(*args, **kwargs)
        # set shifts
        if self.dataset_inst.is_data:
            raise Exception("GetScaleFactorWeights task should only run for MC.")

        if self.normalize_cerrs:
            self.shifts = {"{}_{}".format(shift, direction) for shift, direction in itertools.product(
                ["c_stats1", "c_stats2"], ["up", "down"])}
        else:
            self.shifts = {"nominal"} | {"jes{}_{}".format(shift, direction) for shift, direction in itertools.product(
                jes_sources, ["up", "down"])} | {"{}_{}".format(shift, direction) for shift, direction in itertools.product(
                ["lf", "hf", "lf_stats1", "lf_stats2", "hf_stats1", "hf_stats2"], ["up", "down"])}

    def workflow_requires(self):
        from analysis.tasks.measurement import FitScaleFactors

        reqs = super(GetScaleFactorWeights, self).workflow_requires()

        if not self.cancel_jobs and not self.cleanup_jobs:
            reqs["meta"] = MergeMetaData.req(self, version=self.get_version(MergeMetaData),
                _prefer_cli=["version"])
            reqs["pu"] = CalculatePileupWeights.req(self)
            if not self.pilot:
                reqs["tree"] = MergeTrees.req(self, cascade_tree=-1,
                    version=self.get_version(MergeTrees), _prefer_cli=["version"])

            reqs["sf"] = {shift: FitScaleFactors.req(self, iteration=self.iteration,
                shift=shift, version=self.get_version(FitScaleFactors), _prefer_cli=["version"])
                for shift in self.shifts}

        return reqs

    def requires(self):
        from analysis.tasks.measurement import FitScaleFactors

        reqs = {
            "tree": MergeTrees.req(self, cascade_tree=self.branch, branch=0,
                version=self.get_version(MergeTrees), _prefer_cli=["version", "workflow"]),
            "meta": MergeMetaData.req(self, version=self.get_version(MergeMetaData),
                _prefer_cli=["version"]),
        }
        reqs["pu"] = CalculatePileupWeights.req(self)
        reqs["sf"] = {shift: FitScaleFactors.req(self, iteration=self.iteration,
            shift=shift, version=self.get_version(FitScaleFactors), _prefer_cli=["version"])
            for shift in self.shifts}
        return reqs

    def store_parts(self):
        c_err_part = "c_errors" if self.normalize_cerrs else "b_and_udsg"
        return super(GetScaleFactorWeights, self).store_parts() + (self.iteration,) + (c_err_part,)

    def output(self):
        return self.wlcg_target("stats_{}.json".format(self.branch))

    def get_jec_identifier(self, shift):
        if shift.startswith("jes"):
            return "_" + shift
        else:
            return ""

    def get_scale_factors(self, inp, shift):
        with inp.load() as sfs:
            sf_hists = {}
            for category in sfs.GetListOfKeys():
                category_dir = sfs.Get(category.GetName())
                hist = category_dir.Get("sf")
                # decouple from open file
                hist.SetDirectory(0)

                sf_hists[category.GetName()] = hist

        btag_var = self.config_inst.get_aux("btagger")["variable"]
        identifier = self.get_jec_identifier(shift)

        def get_value(entry):
            scale_factors = []
            for jet_idx in range(1, 5):
                jet_pt = getattr(entry, "jet{}_pt{}".format(jet_idx, identifier))[0]
                jet_eta = getattr(entry, "jet{}_eta{}".format(jet_idx, identifier))[0]
                jet_flavor = getattr(entry, "jet{}_flavor{}".format(jet_idx, identifier))[0]
                jet_btag = getattr(entry, "jet{}_{}{}".format(jet_idx, btag_var, identifier))[0]

                # stop when number of jets is exceeded
                if jet_flavor < -999.:
                    break

                # find category in which the scale factor of the jet was computed to get correct histogram
                if abs(jet_flavor) == 5:
                    region = "hf"
                elif abs(jet_flavor) == 4:
                    region = "c"
                else:
                    region = "lf"

                if region == "c" and not self.normalize_cerrs:
                    continue
                elif region != "c" and self.normalize_cerrs:
                    continue

                category = get_category(jet_pt, abs(jet_eta), region, phase_space="measure")

                # get scale factor
                sf_hist = sf_hists[category.name]
                bin_idx = sf_hist.FindBin(jet_btag)
                scale_factor = sf_hist.GetBinContent(bin_idx)

                scale_factors.append((category, scale_factor))

            return scale_factors

        return get_value

    @law.decorator.notify
    def run(self):
        import ROOT

        inp = self.input()
        outp = self.output()
        outp.parent.touch(0o0770)

        # get child categories
        categories = []
        channels = [self.config_inst.get_aux("dataset_channels")[self.dataset_inst]] \
            if self.dataset_inst.is_data else self.config_inst.channels.values()
        for channel in channels:
            for category, _, children in channel.walk_categories():
                if not children:
                    categories.append((channel, category))
        categories = list(set(categories))

        # get processes
        if len(self.dataset_inst.processes) != 1:
            raise NotImplementedError("only datasets with exactly one linked process can be"
                " handled, got {}".format(len(self.dataset_inst.processes)))
        processes = list(self.dataset_inst.processes.values())
        process = processes[0]

        # prepare dict for outputs
        # shift -> category -> sum weights/ sum weighted sfs
        output_data = {shift: defaultdict(lambda: defaultdict(float)) for shift in self.shifts}

        # open the input file and get the tree
        with inp["tree"].load("READ", cache=False) as input_file:
            tree = input_file.Get("tree")
            self.publish_message("{} events in tree".format(tree.GetEntries()))

            # identifier for jec shifted variables
            for shift in self.shifts:
                jec_identifier = self.get_jec_identifier(shift)

                # pt aliases for jets
                for obj in ["jet1", "jet2", "jet3", "jet4"]:
                    tree.SetAlias("{0}_pt{1}".format(obj, jec_identifier),
                        "({0}_px{1}**2 + {0}_py{1}**2)**0.5".format(obj, jec_identifier))
                # b-tagging alias
                btag_var = self.config_inst.get_aux("btagger")["variable"]
                for obj in ["jet1", "jet2", "jet3", "jet4"]:
                    variable = self.config_inst.get_variable("{0}_{1}".format(obj, btag_var))
                    tree.SetAlias(variable.name + jec_identifier, variable.expression.format(
                        **{"jec_identifier": jec_identifier}))
            # pt aliases for leptons
            for obj in ["lep1", "lep2"]:
                tree.SetAlias("{0}_pt".format(obj),
                    "({0}_px**2 + {0}_py**2)**0.5".format(obj))

            scale_factor_getters = {}
            for shift in self.shifts:
                scale_factor_getters[shift] = self.get_scale_factors(inp["sf"][shift]["sf"], shift)

            # get info to scale event weight to lumi
            x_sec = process.get_xsec(self.config_inst.campaign.ecm).nominal
            sum_weights = inp["meta"].load()["event_weights"]["sum"]
            lumi_factor = x_sec / sum_weights

            input_file.cd()
            with TreeExtender(tree) as te:
                # unpack all branches
                te.unpack_branch("*")
                for entry in te:
                    # get event weight
                    gen_weight = entry.gen_weight[0]
                    channel_id = entry.channel[0]
                    channel = self.config_inst.get_channel(channel_id)
                    lumi = self.config_inst.get_aux("lumi")[channel]

                    evt_weight = gen_weight * lumi * lumi_factor

                    for shift in self.shifts:
                        # event has to pass base selection
                        jec_identifier = self.get_jec_identifier(shift)
                        if getattr(entry, "jetmet_pass{}".format(jec_identifier))[0] != 1:
                            continue

                        # calculate per-jet b-tagging weights
                        scale_factors = scale_factor_getters[shift](entry)
                        # save sum for latter normalization
                        for category, sf_value in scale_factors:
                            output_data[shift][category.name]["sum_sf"] += sf_value * evt_weight
                            output_data[shift][category.name]["sum_weights"] += evt_weight

        # save outputs
        self.output().dump(output_data, formatter="json", indent=4)


class GetScaleFactorWeightsWrapper(WrapperTask):

    wrapped_task = GetScaleFactorWeights


class MergeScaleFactorWeights(AnalysisTask):

    iteration = GetScaleFactorWeights.iteration
    normalize_cerrs = GetScaleFactorWeights.normalize_cerrs

    def requires(self):
        return {dataset.name: GetScaleFactorWeights.req(self, dataset=dataset.name,
            version=self.get_version(MergeScaleFactorWeights), _prefer_cli=["version"])
            for dataset in self.config_inst.datasets if not dataset.is_data}

    def output(self):
        return self.wlcg_target("stats.json")

    def store_parts(self):
        c_err_part = "c_errors" if self.normalize_cerrs else "b_and_udsg"
        return super(MergeScaleFactorWeights, self).store_parts() + (self.iteration,) + (c_err_part,)

    @law.decorator.notify
    def run(self):
        # shift -> category -> values
        stats = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

        def load(inp):
            inp_dict = inp.load()
            for shift, shift_dict in inp_dict.items():
                for category, cat_dict in shift_dict.items():
                    for key, value in cat_dict.items():
                        stats[shift][category][key] += value

        def callback(i):
            self.publish_message("loading meta data {} / {}".format(i + 1, len(input_files)))
            self.publish_progress(100. * (i + 1) / len(input_files))

        input_files = []
        for dataset, inputs in self.input().items():
            input_files.extend(inputs["collection"].targets.values())

        law.util.map_verbose(load, input_files, every=25, callback=callback)

        output_dict = defaultdict(dict)
        for shift, shift_dict in stats.items():
            for category, cat_dict in shift_dict.items():
                output_dict[shift][category] = cat_dict["sum_weights"] / cat_dict["sum_sf"]

        self.output().dump(output_dict, formatter="json", indent=4)
