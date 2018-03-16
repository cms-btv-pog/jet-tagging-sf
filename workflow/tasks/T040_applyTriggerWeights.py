#!/usr/bin/env python

from workflow.base_tasks.base import Task
from workflow.enums import LeptonType, FinalState
from workflow.legacy_tasks.csvSF_TreeReader import CSVSF_TreeReader

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


class ApplyTriggerWeights(Task, law.WrapperTask):
    inclusiveSelection = True
    useTriggerWeights = True

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
                if leptonType == LeptonType.MuonEG and not hf:
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
