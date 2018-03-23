#!/usr/bin/env python
from workflow.enums import Tagger
from workflow.legacy_tasks.rootScriptRunner import RootScriptRunner
from workflow.legacy_tasks.hadd import HAdd

import luigi
import os


class Fit_csvSF(RootScriptRunner):
    tagger = luigi.EnumParameter(
        enum=Tagger,
        default=Tagger.csv,
    )
    isHF = luigi.BoolParameter()
    versionNum = luigi.IntParameter(
        default=0,
    )
    JES = luigi.Parameter(
        default='',
    )

    @property
    def list_scripts(self):
        flavor_str = 'LF'
        if self.isHF:
            flavor_str = 'HF'
        result = ['root_legacy/csvReweightingRun2/csvTreeMaker/'
                  'macros/head13TeV.C']
        result.append(
            "root_legacy/csvReweightingRun2/csvTreeMaker/"
            "macros/fit%s_csvSF_13TeV.C"
            "'(%d, \"%s\", %d, \"%s\", \"%s/\")'"
            % (flavor_str, self.tagger.value, self.input().path,
               self.versionNum, self.JES,
               os.path.dirname(self.output()['root'].path)))
        return result

    @property
    def hist_names(self):
        result = []
        maxPt = 4
        maxEta = 3
        if self.isHF:
            maxPt = 5
            maxEta = 1
        for iPt in range(maxPt):
            for iEta in range(maxEta):
                result.append('csv_ratio_Pt%d_Eta%d' % (iPt, iEta))
        result += [
            'csv_ratio_all',
            # 'csv_ratio_cumulative_all',
        ]
        return result

    def requires(self):
        return HAdd(tagger=self.tagger, versionNum=self.versionNum,
                    JES=self.JES, isHF=self.isHF, version=self.version)

    def output(self):
        flavor_str = 'lf'
        if self.isHF:
            flavor_str = 'hf'
        result = dict()
        dir_name = '%s_v%d%s' % (self.tagger.name, self.versionNum, self.JES)
        result['png'] = []
        for hist_name in self.hist_names:
            for end in ['only', 'JES', 'Stats1', 'Stats2', 'All']:  # 'HF', 'LF'
                result['png'].append(self.local_target(
                    dir_name, 'png', '%sSF_%s_fit_%s.png'
                                     % (flavor_str, hist_name, end)))
            # for end in ['Stats1', 'Stats2']:
            #     result['png'].append(self.local_target(
            #         dir_name, 'png', '%sSF_ratio_%s_fit_%s.png'
            #                          % (flavor_str, hist_name, end)))
        result['root'] = self.local_target(dir_name,
                                           '%s_rwt_fit_%s_v%d%s.root'
                                           % (self.tagger.name, flavor_str,
                                              self.versionNum, self.JES))

        return result
