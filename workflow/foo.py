#!/usr/bin/env python

import os
import subprocess

import luigi
import law

from workflow.base import Task


class FooTask(Task):

    def output(self):
        return self.wlcg_target("foo.json")

    def run(self):
        output = self.output()
        output.parent.touch()
        output.dump({"key": "value"}, indent=4)


# class CalledProcessInfoError(Exception):
#     def __init__(self, message):
#         self.message = message


# class RootScriptRunner(law.Task):
#     batch_mode = luigi.BoolParameter(
#         default=False,
#         significant=False,
#         description='Root Cling flag for batch_mode.',
#     )
#     quit = luigi.BoolParameter(
#         default=False,
#         significant=False,
#         description='Root Cling flag to quit root after execution.',
#     )

#     def __init__(self, list_scripts, *args, **kwargs):
#         super(RootScriptRunner, self).__init__(*args, **kwargs)
#         self.list_scripts = list_scripts

#     def run(self):
#         assert not self.complete()
#         self.clean()

#         exec_comm = 'root'
#         if self.batch_mode:
#             exec_comm += ' -b'
#         if self.quit:
#             exec_comm += ' -q'
#         try:
#             subprocess.check_output([exec_comm] + self.list_scripts,
#                                     stderr=subprocess.STDOUT)
#             assert self.complete()
#         except subprocess.CalledProcessError as e:
#             raise CalledProcessInfoError(
#                 'Command \'%s\' returned non-zero exit status %d.\n%s'
#                 % (e.cmd, e.returncode, e.output.decode('utf-8')))

#     def clean(self):
#         for out in self.output():
#             out.remove()
