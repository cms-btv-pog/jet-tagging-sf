#!/usr/bin/env python

import utils.email
from workflow.base_tasks.base import Task
from workflow.enums import LeptonType, Tagger
from workflow.legacy_tasks.rootScriptRunner import RootScriptRunner
from workflow.tasks.T020_triggerWeights import TriggerWeights

import luigi
import os

# Make plots with checkVar:
# -This produces plots from the root files made from reweighting
#     -Note that the weights are calculated at the end of each iteration,
#      but the plots are made at the beginning of the iteration
#     -So if you want plots with v0 SFs, you need to run the first part of
#      v1 and run checkVar on that
# -Change output directory prefix (L40)
# -Make dirstr point to the directory of files to run on (L103)
# -Change input file suffix for MC and data (e.g. v0 or v3, maybe _ee!)
#  (L106 & L116)
# -Only run on desired data samples (L117 & L172)
# -If running over one of 3 datasets, comment out adding together and just
#  use the desired one (L176-L179)
# -Run with:
#     root -l checkVar.C'(cMVA(0)/csv(1), LF(0)/HF(1), "<suffix>")'
#     -You probably are just doing csv
#     -Run both LF and HF (they can both go to the same output dir)
#         -MuonEG only has HF though!
#     -suffix example: _eeTest


class TriggerWeightsPlots(Task):
    def requires(self):
        requs = dict()
        for leptonType in LeptonType:
            if leptonType == LeptonType.All:
                continue
            if leptonType.name not in requs:
                requs[leptonType.name] = dict()
            for hf, hf_str in zip([True, False], ['hf', 'lf']):
                if leptonType == LeptonType.MuonEG and not hf:
                    continue
                requs[leptonType.name][hf_str] = CheckVarTriggerWeights(
                    isHF=hf,
                    leptonType=leptonType,
                    version=self.version
                )
        return requs

    def output(self):
        return self.local_target('dataMCratio.json')

    def run(self):
        output = self.output()
        output.parent.touch()
        ratios = dict(self.requires())
        for leptonType in LeptonType:
            if leptonType == LeptonType.All:
                continue
            for hf, hf_str in zip([True, False], ['hf', 'lf']):
                if leptonType == LeptonType.MuonEG and not hf:
                    continue
                with open(self.input()[leptonType.name][hf_str]['ratio'].path,
                          'r') as f:
                    ratios[leptonType.name][hf_str] = float(f.readline())
        output.dump(ratios)
        utils.email.notify('TriggerWeightsPlots done')


class CheckVar(RootScriptRunner):
    tagger = luigi.EnumParameter(
        enum=Tagger,
        default=Tagger.csv,
    )
    isHF = luigi.BoolParameter(
        default=True,
    )
    versionNum = luigi.IntParameter(
        default=0,
    )
    leptonType = luigi.EnumParameter(enum=LeptonType)

    @property
    def list_scripts(self):
        result = [
            "root_legacy/csvReweightingRun2/csvTreeMaker/"
            "CSVHistoFiles/checkVar.C"
            "'(%d, %d, \"%s/\", \"%s/\", \"_v%d_histo_All.root\", \"%s\")'"
            % (self.tagger.value, self.isHF,
               os.path.dirname(self.output()['png'][0].path),
               os.path.dirname(self.input()[0].path), self.versionNum,
               self.leptonType.name)
        ]
        return result


class CheckVarTriggerWeights(CheckVar):
    var_list = [
        'probe_jet_csv',
        'nJets',
        'nTags',
        'all_jet_csv',
        'probe_jet_pt',
        'all_jet_pt',
        'first_jet_pt',
        'first_jet_eta',
        'first_jet_csv',
        'first_jet_flavour',
        'second_jet_pt',
        'second_jet_eta',
        'second_jet_csv',
        'second_jet_flavour',
        'met_pt',
        'mht_pt',
        'dr_leplep',
        'mass_leplep',
        '1stlep_pt',
        '2ndlep_pt',
        'lep_pt',
    ]

    def requires(self):
        return TriggerWeights(version=self.version, leptonType=self.leptonType)

    def output(self):
        flavor_str = 'lf'
        if self.isHF:
            flavor_str = 'hf'

        result = []
        for varName in self.var_list:
            fname = '%s_%s_%s.png' % (self.tagger.name, varName, flavor_str)
            result.append(self.local_target(self.leptonType.name, fname))
        result = {'png': result}
        result['ratio'] = self.local_target(self.leptonType.name,
                                            flavor_str + '_dataMCratio.txt')
        return result
