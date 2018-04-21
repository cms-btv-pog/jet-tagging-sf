#!/usr/bin/env python

import utils.email
from workflow.base_tasks.base import Task
from workflow.enums import Tagger
from workflow.legacy_tasks.fit_csvSF import Fit_csvSF

import law

# For CSV:
# (now inclusiveFlag=0) -Change selection (L697) (numjet <2 -> !=2)
# (now inclusiveFlag=0) -Uncomment if(!tpj) block (L897) (do ntag cuts)
# (now inclusiveFlag=0) -Change suffix on root file (histo.root is mandatory)
#                        (L266)
# -Do iterations 0,1,2, and set iter < 3 (runCSV)
# -Run over all data (runCSV at top)
# -The way this works is:
#     -Takes data/mc as input
#     -Applies (or doesn't) selection cuts
#     -Applies previous iteration scale factors (so v0 has none)
#     -Calculates scale factors (remember, these aren't applied until next
#      iteration)


class CSV_SFs(Task, law.WrapperTask):
    def requires(self):
        requs = []
        requs.append(Fit_csvSF(tagger=Tagger.csv, isHF=True, versionNum=2,
                               JES='', version=self.version))
        requs.append(Fit_csvSF(tagger=Tagger.csv, isHF=False, versionNum=2,
                               JES='', version=self.version))
        return requs

    def output(self):
        return self.input()

    def run(self):
        utils.email.notify('CSV_SFs done')
