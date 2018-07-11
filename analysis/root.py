# -*- coding: utf-8 -*-

import ROOT


class ROOTPlot(object):

    tcolors = [1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 23, 28, 30, 42, 48, 50, 95]

    def __init__(self, *args, **kwargs):
        self.canvas = ROOT.TCanvas(*args, **kwargs)

        self.legend = ROOT.TLegend(0.5, 0.7, 0.98, 0.9)
        self.legend.SetNColumns(2)

        self.objects = []

    def draw(self, objs, stacked=False, invis=False, options=[]):
        if not isinstance(options, (list, tuple)):
            options = [options]

        if invis:
            for obj in objs.values():
                invis_obj = obj.Clone()
                self.objects.append(invis_obj)
                invis_obj.SetLineColor(0)
                invis_obj.Draw(" ".join(options))
            return None

        if stacked:
            stack = ROOT.THStack("stack", "")
            for i, (key, obj) in enumerate(objs.items()):
                obj.SetFillColor(self.tcolors[i])
                obj.SetLineColor(self.tcolors[i])
                obj.SetFillStyle(1001)

                self.legend.AddEntry(obj, key, "l")

                stack.Add(obj)
                self.objects.append(obj)
            draw_objs = [stack]
            options.append("HIST")
        else:
            for key, obj in objs.items():
                self.legend.AddEntry(obj, key, "l")
            draw_objs = objs.values()

        for obj in draw_objs:
            self.objects.append(obj)
            obj.Draw(" ".join(options))

    def save(self, path):
        self.legend.Draw()
        self.canvas.Update()
        self.canvas.SaveAs(path)
