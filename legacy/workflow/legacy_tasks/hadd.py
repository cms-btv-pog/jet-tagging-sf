#!/usr/bin/env python
from workflow.enums import Tagger, FinalState, LeptonType
from workflow.legacy_tasks.shellRunner import ShellRunner
from workflow.legacy_tasks.csvSF_TreeReader import CSVSF_TreeReader

import luigi


class HAdd(ShellRunner):
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
    JES = luigi.Parameter(
        default='',
    )

    @property
    def exec_comm(self):
        exec_comm = 'hadd -f %s' % self.output().path
        for inp in self.input():
            exec_comm += ' %s' % inp.path
        return exec_comm

    def requires(self):
        result = []
        for leptonType in LeptonType:
            if leptonType == LeptonType.All:
                continue
            result.append(CSVSF_TreeReader(
                useTriggerWeights=True,
                tagger=self.tagger,
                isHF=self.isHF,
                versionNum=0,
                JES=self.JES,
                insample_ID=leptonType.value,
                leptonType=LeptonType.All,
                inclusiveSelection=False,
                version=self.version))
        for finalState in FinalState:
            result.append(CSVSF_TreeReader(
                useTriggerWeights=True,
                tagger=self.tagger,
                isHF=self.isHF,
                versionNum=self.versionNum,
                JES=self.JES,
                insample_ID=finalState.value,
                leptonType=LeptonType.All,
                inclusiveSelection=False,
                version=self.version))
        return result

    def output(self):
        flavor_str = 'lf'
        if self.isHF:
            flavor_str = 'hf'
        return self.local_target('%s_rwt_%s_all_v%d%s.root'
                                 % (self.tagger.name, flavor_str,
                                    self.versionNum, self.JES))
