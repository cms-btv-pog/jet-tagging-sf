#!/usr/bin/env python

from workflow.base_tasks.base import Task
from workflow.base_tasks.rootScriptRunner import RootScriptRunner
from workflow.enums import LeptonType, RunEra, FinalState, Tagger
from workflow.tasks.T010_loadFiles import Load_Data, Load_MC

import law
import luigi

# For making new trigger weights:
# csvSF_treeReader_13TeV.C:
# -Comment out trigger weights (L628 block)
# -Run over each lep selection separately
#     -Change L839 to the current sample (a,b,c)
# -Change output filename (L277) (e.g. add "_ee")
# -Turn on inclusiveSelection flag
#     -The following is controlled by the inclusiveSelection flag now:
#         -Change selection to inclusive (L718) (numjet !=2 -> <2  )
#         -Comment if(!tpj) block (L922) (don't do ntag cuts)
# runCSVSF.csh:
# -Only do 0 iteration (L23)
# -Set iter > 3 (L34)
#     -This is so that you only run on data, not MC
# -(Un)Comment input datasets (at top) to match csvSF
#
# -Then run runCSVSF.csh
#     -(Once for each lep selection)
#     -This takes awhile (20ish minutes)
#     -Output is in CSVHistoFiles directory


class TriggerWeights(Task, law.WrapperTask):
    inclusiveSelection = True
    useTriggerWeights = False
    leptonType = luigi.EnumParameter(enum=LeptonType)

    def requires(self):
        requs = []
        for hf in [True, False]:
            for finalState in FinalState:
                requs.append(CSVSF_TreeReader(
                    useTriggerWeights=self.useTriggerWeights,
                    isHF=hf,
                    insample_ID=finalState.value,
                    leptonType=self.leptonType,
                    inclusiveSelection=self.inclusiveSelection,
                    version=self.version
                ))
            if self.leptonType == LeptonType.MuonEG and not hf:
                continue
            requs.append(CSVSF_TreeReader(
                useTriggerWeights=self.useTriggerWeights,
                isHF=hf,
                insample_ID=self.leptonType.value,
                leptonType=self.leptonType,
                inclusiveSelection=self.inclusiveSelection,
                version=self.version
            ))
        return requs

    def output(self):
        return self.input()


class CSVSF_TreeReader(RootScriptRunner):
    useTriggerWeights = luigi.BoolParameter()
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
    insample_ID = luigi.IntParameter(
        default=1,
    )
    leptonType = luigi.EnumParameter(enum=LeptonType)
    inclusiveSelection = luigi.BoolParameter()

    @property
    def list_scripts(self):
        result = ['root_legacy/csvReweightingRun2/csvTreeMaker/'
                  'macros/head13TeV.C']
        result.append(
            "root_legacy/csvReweightingRun2/csvTreeMaker/"
            "macros/csvSF_treeReader_13TeV.C"
            "'(\"%s\", \"%s\", %d, %d, %d, %d, %d, \"%s\", %d)'"
            % (self.leptonType.name, self.output().path, self.useTriggerWeights,
               self.inclusiveSelection, self.tagger.value, self.isHF,
               self.versionNum, self.JES, self.insample_ID))
        return result

    def requires(self):
        if self.insample_ID in [i.value for i in LeptonType]:
            requs = []
            for runEra in RunEra:
                requs.append(Load_Data(leptonType=LeptonType(self.insample_ID),
                                       runEra=runEra))
            return requs
        elif self.insample_ID in [i.value for i in FinalState]:
            return Load_MC(finalState=FinalState(self.insample_ID))
        else:
            raise ValueError("%d is not a valid sample number!"
                             % self.insample_ID)

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
        return self.local_target(self.leptonType.name, fname)
