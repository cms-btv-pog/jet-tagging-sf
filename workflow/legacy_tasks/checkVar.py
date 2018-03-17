#!/usr/bin/env python
from workflow.legacy_tasks.rootScriptRunner import RootScriptRunner
from workflow.enums import Tagger, LeptonType

import luigi
import os


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
