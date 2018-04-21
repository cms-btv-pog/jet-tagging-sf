#!/usr/bin/env python
from workflow.legacy_tasks.shellRunner import ShellRunner

import luigi


class RootScriptRunner(ShellRunner):
    batch_mode = luigi.BoolParameter(
        default=True,
        significant=False,
        description='Root Cling flag for batch_mode.',
    )
    quit = luigi.BoolParameter(
        default=True,
        significant=False,
        description='Root Cling flag to quit root after execution.',
    )

    @property
    def exec_comm(self):
        exec_comm = 'root'
        if self.batch_mode:
            exec_comm += ' -b'
        if self.quit:
            exec_comm += ' -q'
        exec_comm += ' ' + ' '.join(self.list_scripts)
        return exec_comm
