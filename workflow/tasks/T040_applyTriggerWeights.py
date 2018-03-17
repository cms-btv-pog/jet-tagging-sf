#!/usr/bin/env python

import utils.email
from workflow.base_tasks.base import Task
from workflow.enums import LeptonType, FinalState
from workflow.legacy_tasks.csvSF_TreeReader import CSVSF_TreeReader
from workflow.legacy_tasks.checkVar import CheckVar

import law

# After getting all the trigger weights, check that things look good:
# csvSF...:
# -Uncomment trigger weights (L628)
# -Change .root file suffix (L277)
# -Use all samples, so comment out override in L839
# runCSV...:
# -Uncomment all 3 data sets (L11)
# checkVar:
# -Change .root file suffix (L104 & L114)
# -Uncomment all fileData (L115 & L170)
# -Make it uses all 3 added together (L174-L178)
#
# -Then make lf & hf plots and post
#     root -l checkVar.C'(cMVA(0)/csv(1), LF(0)/HF(1), "<suffix>")'


class ApplyTriggerWeights(Task, law.WrapperTask):
    inclusiveSelection = True
    useTriggerWeights = True    # TODO right now this just uses the already hard coded trigger weights

    def requires(self):
        requs = []
        for hf in [True, False]:
            for finalState in FinalState:
                requs.append(CSVSF_TreeReader(
                    useTriggerWeights=self.useTriggerWeights,
                    isHF=hf,
                    insample_ID=finalState.value,
                    leptonType=LeptonType.All,
                    inclusiveSelection=self.inclusiveSelection,
                    version=self.version
                ))
            for leptonType in LeptonType:
                if leptonType == LeptonType.All:
                    continue
                requs.append(CSVSF_TreeReader(
                    useTriggerWeights=self.useTriggerWeights,
                    isHF=hf,
                    insample_ID=leptonType.value,
                    leptonType=LeptonType.All,
                    inclusiveSelection=self.inclusiveSelection,
                    version=self.version
                ))
        return requs

    def output(self):
        return self.input()


class AppliedTriggerWeightsPlots(Task):
    def requires(self):
        requs = dict()
        for hf, hf_str in zip([True, False], ['hf', 'lf']):
            requs[hf_str] = CheckVarAppliedTriggerWeights(
                isHF=hf,
                version=self.version
            )
        return requs

    def output(self):
        return self.local_target('dataMCratio.json')

    def run(self):
        output = self.output()
        output.parent.touch()
        ratios = dict()
        for hf, hf_str in zip([True, False], ['hf', 'lf']):
            with open(self.input()[hf_str]['ratio'].path,
                      'r') as f:
                ratios[hf_str] = float(f.readline())
        output.dump(ratios)
        utils.email.notify('AppliedTriggerWeightsPlots done')


class CheckVarAppliedTriggerWeights(CheckVar):
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
    leptonType = LeptonType.All

    def requires(self):
        return ApplyTriggerWeights(version=self.version)

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
