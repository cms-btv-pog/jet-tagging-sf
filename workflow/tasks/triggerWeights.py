#!/usr/bin/env python

from workflow.base_tasks.base import Task
from workflow.base_tasks.rootScriptRunner import RootScriptRunner

import enum
import law
import luigi
import os

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


class Tagger(enum.IntEnum):
    cMVA = 0
    csv = 1


class LeptonType(enum.IntEnum):
    DoubleEG = -100
    DoubleMuon = -200
    MuonEG = -300


class RunEra(enum.Enum):
    B = 1
    C = 2
    D = 3
    E = 4
    F = 5


class FinalState(enum.IntEnum):
    zjets = 2300
    lowMasszjets = 2310
    # WJetsToLNu = 2400
    ttjets = 2500
    # TToLepton_s = 2510
    # TBarToLepton_s = 2511
    # TToLeptons_t = 2512
    # TBarToLeptons_t = 2513
    singletW = 2514
    singletbarW = 2515
    # TTZJets = 2523
    # TTWJets = 2524
    WW = 2600


class Load_Data(luigi.ExternalTask, Task):
    version = None
    leptonType = luigi.EnumParameter(enum=LeptonType)
    runEra = luigi.EnumParameter(enum=RunEra)
    date = luigi.Parameter(
        default='2018-02-10',
    )

    def output(self):
        return law.LocalFileTarget(os.path.join(
            self.local_data_root,
            'in',
            self.date,
            '%s_Run2017%s_17Nov2017_%s.root'
            % (self.leptonType.name, self.runEra.name, self.date)
        ))


class Load_MC(luigi.ExternalTask, Task):
    version = None
    finalState = luigi.EnumParameter(enum=FinalState)
    date = luigi.Parameter(
        default='2018-02-10',
    )

    def output(self):
        return law.LocalFileTarget(os.path.join(
            self.local_data_root,
            'in',
            self.date,
            '%s_%s.root' % (self.finalState.name, self.date)
        ))


class TriggerWeights(Task, law.WrapperTask):
    inclusiveSelection = True
    useTriggerWeights = False

    def requires(self):
        requs = []
        for leptonType in LeptonType:
            for hf in [True, False]:
                if leptonType == LeptonType.MuonEG and not hf:
                    continue
                requs.append(CSVSF_TreeReader(
                    useTriggerWeights=self.useTriggerWeights,
                    isHF=hf,
                    insample_ID=leptonType.value,
                    inclusiveSelection=self.inclusiveSelection,
                    version=self.version
                ))
        for finalState in FinalState:
            for hf in [True, False]:
                requs.append(CSVSF_TreeReader(
                    useTriggerWeights=self.useTriggerWeights,
                    isHF=hf,
                    insample_ID=finalState.value,
                    inclusiveSelection=self.inclusiveSelection,
                    version=self.version
                ))
        return requs


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
    inclusiveSelection = luigi.BoolParameter()

    @property
    def list_scripts(self):
        result = ['root_legacy/csvReweightingRun2/csvTreeMaker/'
                  'macros/head13TeV.C']
        result.append(
            "root_legacy/csvReweightingRun2/csvTreeMaker/"
            "macros/csvSF_treeReader_13TeV.C"
            "'(\"%s\", %d, %d, %d, %d, %d, \"%s\", %d)'"
            % (self.output().path, self.useTriggerWeights,
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
        return self.local_target(fname)
