# -*- coding: utf-8 -*-

import ROOT


class ROOTPlot(object):

    def __init__(self, *args, **kwargs):
        ROOT.PyConfig.IgnoreCommandLineOptions = True
        ROOT.gROOT.SetBatch()

        self.canvas = ROOT.TCanvas(*args, **kwargs)
        self.objects = []

    def draw(self, objs, stacked=False, options=[]):
        if not isinstance(objs, (list, tuple)):
            objs = [objs]
        if not isinstance(options, (list, tuple)):
            options = [options]

        if stacked:
            stack = ROOT.THStack("stack", "")
            for i, obj in enumerate(objs):
                obj.SetFillColor(i + 1)
                obj.SetLineColor(i + 1)
                obj.SetFillStyle(1001)

                stack.Add(obj)
                self.objects.append(obj)
            objs = [stack]
            options.append("HIST")

        for obj in objs:
            self.objects.append(obj)
            obj.Draw(" ".join(options))

    def save(self, path):
        self.canvas.Update()
        self.canvas.SaveAs(path)
