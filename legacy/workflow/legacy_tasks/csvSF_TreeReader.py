#!/usr/bin/env python
from workflow.enums import LeptonType, RunEra, FinalState, Tagger
from workflow.legacy_tasks.rootScriptRunner import RootScriptRunner
from workflow import legacy_tasks as lt
from workflow.tasks.T010_loadFiles import Load_Data, Load_MC

import luigi
import os


class CSVSF_TreeReader(RootScriptRunner):
    useTriggerWeights = luigi.BoolParameter()
    tagger = luigi.EnumParameter(
        enum=Tagger,
        default=Tagger.csv,
    )
    isHF = luigi.BoolParameter()
    versionNum = luigi.IntParameter(
        default=0,
    )
    JES = luigi.Parameter(
        default='',
    )
    insample_ID = luigi.IntParameter(
        default=1,
    )
    leptonType = luigi.EnumParameter(enum=LeptonType)
    inclusiveSelection = luigi.BoolParameter()

    @property
    def list_scripts(self):
        if self.versionNum == 0:
            prevVerDir = ''
        else:
            prevVerDir = os.path.dirname(
                self.input()['prevVer'][0]['root'].path)
        result = ['root_legacy/csvReweightingRun2/csvTreeMaker/'
                  'macros/head13TeV.C']
        result.append(
            "root_legacy/csvReweightingRun2/csvTreeMaker/"
            "macros/csvSF_treeReader_13TeV.C"
            "'(\"%s\", \"%s\", \"%s\", %d, %d, %d, %d, %d, \"%s\", %d)'"
            % (prevVerDir, self.leptonType.name, self.output().path,
               self.useTriggerWeights, self.inclusiveSelection,
               self.tagger.value, self.isHF, self.versionNum, self.JES,
               self.insample_ID))
        return result

    def requires(self):
        if self.insample_ID in map(int, LeptonType):
            requs_in = []
            for runEra in RunEra:
                requs_in.append(Load_Data(
                    leptonType=LeptonType(self.insample_ID),
                    runEra=runEra))
        elif self.insample_ID in [i.value for i in FinalState]:
            requs_in = Load_MC(finalState=FinalState(self.insample_ID))
        else:
            raise ValueError("%d is not a valid sample number!"
                             % self.insample_ID)
        if self.versionNum > 0:
            requs_prevVer = []
            requs_prevVer.append(lt.fit_csvSF.Fit_csvSF(
                tagger=self.tagger, isHF=True, version=self.version,
                versionNum=self.versionNum - 1, JES=self.JES))
            requs_prevVer.append(lt.fit_csvSF.Fit_csvSF(
                tagger=self.tagger, isHF=False, version=self.version,
                versionNum=self.versionNum - 1, JES=self.JES))
            return {'prevVer': requs_prevVer,
                    'in': requs_in}
        else:
            return {'in': requs_in}

    def output(self):
        if self.insample_ID in map(int, LeptonType):
            sample_name = LeptonType(self.insample_ID)
        elif self.insample_ID in map(int, FinalState):
            sample_name = FinalState(self.insample_ID)
        else:
            raise ValueError("%d is not a valid sample number!"
                             % self.insample_ID)
        flavor_str = 'lf'
        if self.isHF:
            flavor_str = 'hf'
        fname = '%s_rwt_%s_%s_v%d%s' % (self.tagger.name, flavor_str,
                                        sample_name.name, self.versionNum,
                                        self.JES)
        if self.inclusiveSelection:
            fname += '_histo_All.root'
        else:
            fname += '_histo.root'
        trigger_str = ''
        if self.useTriggerWeights:
            trigger_str = '_withTrigWgts'
        return self.local_target(self.leptonType.name + trigger_str, fname)
