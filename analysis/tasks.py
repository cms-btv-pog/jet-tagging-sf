# -*- coding: utf-8 -*-

import os
import subprocess
import collections

import law
import luigi

from analysis.base import AnalysisTask, DatasetTask, GridWorkflow


class GetDatasetLFNs(DatasetTask, law.TransferLocalFile):

    source_path = None
    version = None

    def single_output(self):
        h = law.util.create_hash(sorted(self.dataset_inst.keys))
        return self.wlcg_target("lfns_{}.json".format(h))

    def run(self):
        lfns = []
        for key in self.dataset_inst.keys:
            print("get lfns for key {}".format(key))
            cmd = "dasgoclient -query='file dataset={}' -limit=0".format(key)
            code, out, _ = law.util.interruptable_popen(cmd, shell=True, stdout=subprocess.PIPE,
                executable="/bin/bash")
            if code != 0:
                raise Exception("dasgoclient query failed")
            lfns.extend(out.strip().split("\n"))

        tmp = law.LocalFileTarget(is_tmp="json")
        tmp.dump(lfns)
        self.transfer(tmp)


class CreateTuples(DatasetTask, GridWorkflow, law.LocalWorkflow):

    def create_branch_map(self):
        return list(range(self.dataset_inst.n_files))

    def requires(self):
        return {
            "lfns": GetDatasetLFNs.req(self),
        }

    def output(self):
        return self.wlcg_target("tuple_{}.root".format(self.branch_data))

    def run(self):
        pass
