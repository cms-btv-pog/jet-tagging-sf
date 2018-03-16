#!/usr/bin/env python
from workflow.base_tasks.base import Task

import luigi
import subprocess


class RootScriptRunner(Task):
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

    def __init__(self, *args, **kwargs):
        super(RootScriptRunner, self).__init__(*args, **kwargs)

    def run(self):
        assert not self.complete()
        self.clean()

        exec_comm = 'root'
        if self.batch_mode:
            exec_comm += ' -b'
        if self.quit:
            exec_comm += ' -q'
        exec_comm += ' ' + ' '.join(self.list_scripts)
        try:
            print('Running command "%s" ...' % exec_comm)
            subprocess.check_output(exec_comm, shell=True,
                                    stderr=subprocess.STDOUT)
            if not self.complete():
                raise RuntimeError(
                    'The command "%s" did not produce the expected outputs:'
                    '\n%s' % (exec_comm,
                              '\n'.join([i.path for i in
                                         luigi.task.flatten_output(self)])))
        except subprocess.CalledProcessError as e:
            print('Command \'%s\' returned non-zero exit status %d.\n%s'
                  % (e.cmd, e.returncode, e.output.decode('utf-8')))
            raise e
        except OSError as e:
            if e.filename is not None:
                e.strerror = 'No such file or directory %s.' % e.filename
                raise e
            e.strerror = 'No valid command "%s".' % exec_comm
            raise e

    def clean(self):
        for out in luigi.task.flatten_output(self):
            out.remove()
            out.parent.touch()
