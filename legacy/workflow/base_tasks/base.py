#!/usr/bin/env python

import os

import luigi
import law

law.contrib.load("wlcg")


class Task(law.Task):

    version = luigi.Parameter(description="task version, required")

    exclude_db = True

    def __init__(self, *args, **kwargs):
        super(Task, self).__init__(*args, **kwargs)

        # other attributes
        self.local_data_root = os.getenv("JTSF_DATA")

    def store_parts(self):
        return (self.__class__.__name__,)

    def store_parts_opt(self):
        parts = tuple()
        if not law.is_no_param(self.version):
            parts += (self.version,)
        return parts

    def local_store(self):
        parts = (self.local_data_root,) + self.store_parts() \
            + self.store_parts_opt()
        return os.path.join(*parts)

    def local_path(self, *parts):
        return os.path.join(self.local_store(), *[str(part) for part in parts])

    def local_target(self, *parts):
        return law.LocalFileTarget(self.local_path(*parts))

    def wlcg_store(self):
        parts = tuple() + self.store_parts() + self.store_parts_opt()
        return os.path.join(*parts)

    def wlcg_path(self, *parts):
        return os.path.join(self.wlcg_store(), *[str(part) for part in parts])

    def wlcg_target(self, *parts):
        return law.WLCGFileTarget(self.wlcg_path(*parts))
