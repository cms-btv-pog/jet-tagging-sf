# -*- coding: utf-8 -*-


import os
import uuid
import shutil
import subprocess
import collections

import law
import luigi
import six

from analysis.tasks.base import AnalysisTask, DatasetTask
from analysis.util import wget


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


class DownloadSetupFiles(AnalysisTask, law.TransferLocalFile):

    source_path = None
    version = None

    def __init__(self, *args, **kwargs):
        super(DownloadSetupFiles, self).__init__(*args, **kwargs)

        self.source_files = self.get_source_files()

    def single_output(self):
        h = law.util.create_hash(law.util.flatten(self.source_files))
        return self.wlcg_target("{}.tgz".format(h))

    @classmethod
    def create_path_hash(cls, path):
        if not isinstance(path, six.string_types):
            return None
        elif not path.startswith("http") and not path.startswith("/"):
            return None
        else:
            return "{}_{}".format(law.util.create_hash(path), os.path.basename(path))

    def get_source_files(self):
        # prepare JES files
        jes_url = lambda version, level: "https://raw.githubusercontent.com/cms-jet/JECDatabase" \
            "/master/textFiles/{0}/{0}_{1}_AK4PFchs.txt".format(version, level)
        jes_files = collections.defaultdict(lambda: collections.defaultdict(dict))
        for src in ("mc", "data"):
            for _, _, version in self.config_inst.get_aux("jes_version")[src]:
                for level in self.config_inst.get_aux("jes_levels")[src] + ["Uncertainty"]:
                    jes_files[src][version][level] = jes_url(version, level)
        jes_unc_src_file = jes_url(self.config_inst.get_aux("jes_version")["mc"][0][2],
            "UncertaintySources")

        # prepare JER files
        jer_url = lambda version, src: "https://raw.githubusercontent.com/cms-jet/JRDatabase" \
            "/master/textFiles/{0}/{0}_{1}_AK4PFchs.txt".format(version, src)
        jer_files = {}
        for src in ("SF", "PtResolution", "PhiResolution"):
            jer_files[src] = jer_url(self.config_inst.get_aux("jer_version") + "_MC", src)

        return {
            "lumi_file": self.config_inst.get_aux("lumi_file"),
            "normtag_file": self.config_inst.get_aux("normtag_file"),
            "jes_files": jes_files,
            "jes_unc_src_file": jes_unc_src_file,
            "jer_files": jer_files,
        }

    def run(self):
        # create a tmp dir
        tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
        tmp_dir.touch()

        # download all setup files
        def download(src):
            h = self.create_path_hash(src)
            if h is None:
                return
            self.publish_message("download {}".format(src))
            dst = os.path.join(tmp_dir.path, h)
            if src.startswith("http"):
                wget(src, dst)
            else:
                shutil.copy2(src, dst)

        law.util.map_struct(download, self.source_files)

        # create a tmp archive
        tmp_arc = law.LocalFileTarget(is_tmp="tgz")
        tmp_arc.dump(tmp_dir)

        # transfer
        self.transfer(tmp_arc)

    def localize(self, **kwargs):
        # load the archive and unpack it into a temporary directory
        tmp_dir = law.LocalDirectoryTarget(is_tmp=True)
        output = self.output()
        if self.replicas >= 1:
            output = output.random_target()
        output.load(tmp_dir, **kwargs)

        def abspath(path):
            h = self.create_path_hash(path)
            return h and os.path.join(tmp_dir.path, h)

        return tmp_dir, law.util.map_struct(abspath, self.source_files)


class CalculateLumi(AnalysisTask):

    channel = luigi.Parameter(default="ee", description="analysis channel")
    period = luigi.Parameter(default=law.NO_STR, description="character denoting the data "
        "taking period, default: empty")
    version = None

    def __init__(self, *args, **kwargs):
        super(CalculateLumi, self).__init__(*args, **kwargs)

        # get the channel instance
        self.channel_inst = self.config_inst.get_channel(self.channel)

        # get the run range if period is set
        self.run_range = None
        if self.period != law.NO_STR:
            self.run_range = self.config_inst.get_aux("run_ranges")[self.period]

        self.has_run = False

    def complete(self):
        return self.has_run

    def requires(self):
        return DownloadSetupFiles.req(self)

    @law.decorator.notify
    def run(self):
        if not law.util.check_bool_flag(os.getenv("JTSF_ON_LXPLUS")):
            raise Exception("{} must run on lxplus".format(self.__class__.__name__))

        triggers = self.config_inst.get_aux("triggers")[self.channel_inst]
        setup_files_dir, setup_files = self.requires().localize()
        uid = str(uuid.uuid4())

        # a tmp dir
        tmp = law.LocalDirectoryTarget(is_tmp=True)
        tmp.touch()

        # build the command
        triggers_str = " ".join(triggers)
        begin_end = "--begin {} --end {}".format(*self.run_range) if self.run_range else ""
        cmd = """
            export PATH="$( pwd )/bin:/afs/cern.ch/cms/lumi/brilconda/bin:$PATH"
            export PYTHONPATH="$( pwd )/lib/python2.7/site-packages:$PYTHONPATH"
            source activate root
            pip install --prefix . --ignore-installed brilws
            >&2 echo "using brilcalc $( brilcalc --version ) from $( which brilcalc )"
            >&2 echo "lumi file: {lumi_file}"
            >&2 echo "norm file: {normtag_file}"
            >&2 echo "triggers : {triggers}"
            >&2 echo "run range: {begin_end}"
            for HLTPATH in {triggers}; do
                >&2 echo "calculate lumi for trigger path $HLTPATH ..."
                brilcalc lumi \
                    -u /pb \
                    --hltpath "$HLTPATH" \
                    --normtag "{normtag_file}" \
                    -i "{lumi_file}" \
                    -b "STABLE BEAMS" \
                    {begin_end} \
                || exit "$?"
                echo "{uid}"
                >&2 echo "done"
            done
        """.format(lumi_file=setup_files["lumi_file"], normtag_file=setup_files["normtag_file"],
                triggers=triggers_str, begin_end=begin_end, uid=uid)

        # run the command
        code, out, _ = law.util.interruptable_popen(cmd, shell=True, executable="/bin/bash",
            stdout=subprocess.PIPE, cwd=tmp.path)
        if code != 0:
            raise Exception("brilcalc failed")

        # parse the output
        blocks = out.split(uid)[:-1]
        lumi_data = {}
        for trigger, block in zip(triggers, blocks):
            lines = block[:block.find("#Summary")].strip().split("\n")[:-1]

            # traverse backwards until a line does not start with "|"
            # columns: run:fill, time, ncms, hltpath, delivered, recorded
            while lines:
                line = lines.pop().strip()
                if not line.startswith("|"):
                    break

                parts = [p.strip() for p in line.split("|")[1:-1]]
                run = int(parts[0].split(":")[0])
                path = parts[3]
                lumi = float(parts[5])
                lumi_data.setdefault(run, {})[path] = lumi

        # calculate the lumi
        lumi = 0.
        for data in lumi_data.values():
            # data is a dict "hlt path -> lumi" per run
            # multiple elements mean that multiple, OR-connected triggers were active in that run
            # in this case, use the maximum as smaller values result from prescales
            lumi += max(list(data.values()))

        self.publish_message("Integrated luminosity: {} /pb".format(lumi))
