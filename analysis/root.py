# -*- coding: utf-8 -*-

import ROOT

tcolors = {
    "data": ROOT.kBlack,
    "tt_dl": ROOT.kBlue,
    "dy_lep_5To50_Ht70To100": ROOT.kSpring + 1,
    "dy_lep_5To50_Ht100To200": ROOT.kSpring + 10,
    "dy_lep_5To50_Ht200To400": ROOT.kGreen,
    "dy_lep_5To50_Ht400To600": ROOT.kGreen - 9,
    "dy_lep_5To50_Ht600ToInf": ROOT.kTeal,
    "dy_lep_10To50": ROOT.kGreen,
    "dy_lep_50ToInf": ROOT.kTeal - 5,
    "dy_lep_50ToInf_Ht70To100": ROOT.kYellow,
    "dy_lep_50ToInf_Ht100To200": ROOT.kOrange - 3,
    "dy_lep_50ToInf_Ht200To400": ROOT.kOrange + 7,
    "dy_lep_50ToInf_Ht400To600": ROOT.kRed,
    "dy_lep_50ToInf_Ht600To800": ROOT.kPink - 3,
    "dy_lep_50ToInf_Ht800To1200": ROOT.kPink + 7,
    "dy_lep_50ToInf_Ht1200To2500": ROOT.kMagenta,
    "dy_lep_50ToInf_Ht2500ToInf": ROOT.kMagenta - 5,
    "st_tW_t": ROOT.kCyan,
    "st_tW_tbar": ROOT.kCyan - 10,
    "WW_sl": ROOT.kGray,
    "b": ROOT.kRed,
    "c": ROOT.kOrange,
    "udsg": ROOT.kGreen,
}


class ROOTPad(object):
    def __init__(self, *args, **kwargs):
        self.pad = ROOT.TPad(*args, **kwargs)
        self.pad.Draw()

        self.legend = ROOT.TLegend(0.5, 0.7, 0.88, 0.88)
        self.legend.SetNColumns(2)

        self.objects = []

    def draw(self, obj_dict, stacked=False, invis=False, options=[]):
        if not isinstance(options, (list, tuple)):
            options = [options]

        if invis:
            # Draw an invisible object to fix the axis ranges. Required if the first object to be
            # drawn is a THStack
            for obj in obj_dict.values():
                invis_obj = obj.Clone()
                self.objects.append(invis_obj)
                invis_obj.SetLineColor(0)
                invis_obj.SetMarkerColor(0)
                invis_obj.Draw(" ".join(options))
            return None

        if stacked:
            stack = ROOT.THStack("stack", "")
            for key, obj in sorted(obj_dict.items()):
                obj.SetFillColor(tcolors.get(key, 1))
                obj.SetLineColor(tcolors.get(key, 1))
                obj.SetFillStyle(1001)

                self.legend.AddEntry(obj, key, "f")

                stack.Add(obj)
                self.objects.append(obj)
            options.append("HIST")
            draw_objs = [stack]
        else:
            for key, obj in sorted(obj_dict.items()):
                self.legend.AddEntry(obj, key, "l")
                obj.SetLineColor(tcolors.get(key, 1))
            draw_objs = obj_dict.values()

        for obj in draw_objs:
            self.objects.append(obj)
            obj.Draw(" ".join(options))

    def save(self):
        self.pad.cd()
        if hasattr(self, "legend"):
            self.legend.Draw()
        self.pad.Update()

    def close(self):
        self.pad.Close()

    def cd(self):
        self.pad.cd()


class ROOTPlot(object):

    def __init__(self, *args, **kwargs):
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetLegendBorderSize(0)

        self.canvas = ROOT.TCanvas(*args, **kwargs)

        # store all used ROOT object to make sure they are not cleaned up
        self.objects = []

    def create_pads(self, n_pads_x=1, n_pads_y=1, limits_x=None, limits_y=None):
        if limits_x is None:
            limits_x = [idx / float(n_pads_x) for idx in range(n_pads_x + 1)]
        if limits_y is None:
            limits_y = [idx / float(n_pads_y) for idx in range(n_pads_y + 1)]

        self.pads = {}
        for idx_x in xrange(n_pads_x):
            for idx_y in xrange(n_pads_y):
                self.canvas.cd()
                name = "{}_{}_{}".format(hash(self), idx_x, idx_y)
                pad = ROOTPad(name, name, limits_x[idx_x], limits_y[idx_y],
                    limits_x[idx_x + 1], limits_y[idx_y + 1])
                self.pads[(idx_x, idx_y)] = pad
                self.open_pad = pad

    def cd(self, idx_x, idx_y):
        self.canvas.cd()
        self.open_pad = self.pads[(idx_x, idx_y)]
        self.open_pad.cd()

    def draw(self, *args, **kwargs):
        self.open_pad.draw(*args, **kwargs)

    def save(self, path):
        for pad in self.pads.values():
            pad.save()
        self.canvas.Update()
        self.canvas.SaveAs(path)

    def __del__(self):
        for pad in self.pads.values():
            pad.close()
