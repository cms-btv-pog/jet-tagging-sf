#!/usr/bin/env python

from workflow.base_tasks.base import Task
from workflow.enums import LeptonType, FinalState
from workflow.legacy_tasks.csvSF_TreeReader import CSVSF_TreeReader

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
        assert self.leptonType != LeptonType.All
        requs = []
        hfs = [True, False]
        if self.leptonType == LeptonType.MuonEG:
            hfs = [True]
        for hf in hfs:
            for finalState in FinalState:
                requs.append(CSVSF_TreeReader(
                    useTriggerWeights=self.useTriggerWeights,
                    isHF=hf,
                    insample_ID=finalState.value,
                    leptonType=self.leptonType,
                    inclusiveSelection=self.inclusiveSelection,
                    version=self.version
                ))
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
