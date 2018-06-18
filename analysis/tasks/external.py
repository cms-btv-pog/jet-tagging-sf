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
    version = None

    def __init__(self, *args, **kwargs):
        super(CalculateLumi, self).__init__(*args, **kwargs)

        # get the channel instance
        self.channel_inst = self.config_inst.get_channel(self.channel)

        self.has_run = False

    def complete(self):
        return self.has_run

    def requires(self):
        return DownloadSetupFiles.req(self)

    def run(self):
        triggers = self.config_inst.get_aux("triggers")[self.channel_inst]
        setup_files_dir, setup_files = self.requires().localize()
        uid = str(uuid.uuid4())

        # a tmp dir
        tmp = law.LocalDirectoryTarget(is_tmp=True)
        tmp.touch()

        # build the command
        cmd = """
            export PATH="$( pwd )/bin:/afs/cern.ch/cms/lumi/brilconda/bin:$PATH"
            export PYTHONPATH="$( pwd )/lib/python2.7/site-packages:$PYTHONPATH"
            source activate root
            pip install --prefix . --upgrade brilws
            >&2 echo "using brilcalc $( brilcalc --version ) from $( which brilcalc )"
            >&2 echo "lumi file: {lumi_file}"
            >&2 echo "norm file: {normtag_file}"
            >&2 echo "triggers : {triggers}"
            for HLTPATH in {triggers}; do
                >&2 echo "calculate lumi for trigger path $HLTPATH"
                brilcalc lumi \
                    -u /pb \
                    --hltpath "$HLTPATH" \
                    --normtag "{normtag_file}" \
                    -i "{lumi_file}" \
                    -b "STABLE BEAMS" || exit "$?"
                    echo "{uid}"
            done
        """.format(lumi_file=setup_files["lumi_file"], normtag_file=setup_files["normtag_file"],
                triggers=" ".join(triggers[:1]), uid=uid)

        # run the command
        code, out, _ = law.util.interruptable_popen(cmd, shell=True, executable="/bin/bash",
            stdout=None, cwd=tmp.path)
        if code != 0:
            raise Exception("brilcalc failed")

        # print("you might need to enter your lxplus password multiple times")
        # p = Popen(cmd, stdout=PIPE, stderr=sys.stderr, shell=True, executable="/bin/bash")
        # out = p.communicate()[0]
        # if p.returncode != 0:
        #     raise Exception("lumi calculation query failed (%s)" % p.returncode)

        # lumiData = OrderedDict()
        # for trigger, result in zip(triggers, out.strip().split(uid)):
        #     print("trigger: %s" % trigger)
        #     print(result)
        #     idx = result.find("#Summary")

        #     # parse output
        #     lines = result[:idx+1].split("\n")
        #     lineIdxs = [i for i, line in enumerate(lines) if line.startswith("+")]
        #     startLine = lineIdxs[1] + 1
        #     endLine = lineIdxs[-1]
        #     lines = lines[startLine:endLine]
        #     for line in lines:
        #         elems = [elem.strip() for elem in line.split("|")[1:-1]]

        #         run = int(elems[0].split(":")[0])
        #         lumi = float(elems[5])
        #         _trigger = elems[3]
        #         lumiData.setdefault(run, {})[_trigger] = lumi

        # # some summaries
        # ll = LumiList(lumiData)
        # print("\nLuminosities [/pb]:")
        # for trigger in triggers:
        #     print("%s: %.3f" % (trigger, ll.get(hltPath=trigger)))
        # print("Total: %.3f\n" % ll.get())

        # with self.output().localTmp() as tmp:
        #     with tmp.open("w") as f:
        #         json.dump(lumiData, f, indent=4)
