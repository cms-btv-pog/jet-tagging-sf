#!/usr/bin/env python

from workflow.base_tasks.base import Task
from workflow.enums import LeptonType, RunEra, FinalState

import law
import luigi
import os


class Load_Data(Task, law.ExternalTask):
    version = None
    leptonType = luigi.EnumParameter(enum=LeptonType)
    runEra = luigi.EnumParameter(enum=RunEra)
    date = luigi.Parameter(
        default='2018-02-10',
    )

    def output(self):
        assert self.leptonType.value < 0
        return law.LocalFileTarget(os.path.join(
            self.local_data_root,
            'in',
            self.date,
            '%s_Run2017%s_17Nov2017_%s.root'
            % (self.leptonType.name, self.runEra.name, self.date)
        ))


class Load_MC(Task, law.ExternalTask):
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
