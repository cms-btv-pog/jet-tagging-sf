# -*- coding: utf-8 -*-

import ROOT

from array import array

tcolors = {
    "data": ROOT.kBlack,
    "tt_dl": ROOT.kGray,
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

        self.left_margin = self.pad.GetLeftMargin()
        self.top_margin = self.pad.GetTopMargin()
        self.right_margin = self.pad.GetRightMargin()
        self.bottom_margin = self.pad.GetBottomMargin()

        self.objects = []
        self.missing_key = False # Legend contains key not found in tcolors dict
        self.line_colors = [2, 4, 3, 90, 6, 8]
        self.has_drawn_object = False # To automatically set option 'SAME' if needed


    def draw_base_legend(self):
        self.legend = ROOT.TLegend(0.5, 1.2 * self.bottom_margin, 1. - 1.2 * self.right_margin, 0.4)
        self.legend.SetNColumns(2)

    def draw_text(self, text, xpos=None, ypos=None, size=None):
        self.cd()
        if xpos is None:
            xpos = 1. - 1.2 * self.right_margin
        if ypos is None:
            ypos = 2. - 1.2 * self.top_margin
        if size is None:
            size = 0.75 * self.top_margin

        root_text = ROOT.TLatex()
        root_text.SetTextSize(size)
        root_text.SetTextAlign(33)
        root_text.SetTextFont(42)
        root_text.DrawLatex(xpos, ypos, text)
        self.objects.append(root_text)

    def add_object(self, obj):
        # separate objects from their root files so that the plot persists if the file is closed
        if hasattr(obj, "SetDirectory"):
            obj.SetDirectory(0)
        self.objects.append(obj)

    def update_options(self, options, add_same_option):
        if options is None:
            options = []
        if not isinstance(options, (list, tuple)):
            options = [options]

        if add_same_option and self.has_drawn_object and not "SAME" in options:
            options.append("SAME")
        self.has_drawn_object = True

        return options

    def get_line_color(self, line_color):
        if line_color is None:
            line_color = self.line_colors.pop(0)
        return line_color

    def draw(self, obj_dict, stacked=False, invis=False, line_color=1, fill_color=1,
        stack_maximum=None, options=None, add_same_option=True, add_to_legend=True):
        self.cd()
        line_color = self.get_line_color(line_color)

        options = self.update_options(options, add_same_option)

        if invis:
            # Draw an invisible object to fix the axis ranges.
            for obj in obj_dict.values():
                invis_obj = obj.Clone()
                self.add_object(invis_obj)
                invis_obj.SetLineColor(0)
                invis_obj.SetMarkerColor(0)
                invis_obj.Draw(" ".join(options))
            return None

        if stacked:
            from random import randint
            stack = ROOT.THStack("stack_" + str(randint(0, 10**9)), str(randint(0, 10**9)))
            for key, obj in sorted(obj_dict.items(), reverse=True):
                if key not in tcolors:
                    key = "Other"
                    if add_to_legend and not self.missing_key:
                        self.legend.AddEntry(obj, key, "f")
                        self.missing_key = True
                elif add_to_legend:
                    self.legend.AddEntry(obj, key, "f")

                obj.SetFillColor(tcolors.get(key, fill_color))
                obj.SetLineColor(tcolors.get(key, line_color))
                obj.SetFillStyle(1001)

                stack.Add(obj)
                self.add_object(obj)
            options.append("HIST")
            if stack_maximum is not None:
                stack.SetMaximum(stack_maximum)
            draw_objs = [stack]
        else:
            for key, obj in sorted(obj_dict.items()):
                obj.SetLineColor(tcolors.get(key, line_color))
                if add_to_legend:
                    self.legend.AddEntry(obj, key, "l")
            draw_objs = obj_dict.values()

        for obj in draw_objs:
            self.add_object(obj)
            obj.Draw(" ".join(options))

    def draw_as_graph(self, hist, options=None, add_same_option=True):
        self.cd()
        options = self.update_options(options, add_same_option)

        x, y = [], []
        xerr_down, xerr_up = [], []
        yerr_down, yerr_up = [], []

        for i in xrange(1, hist.GetNbinsX() + 1):
            x.append(hist.GetBinCenter(i))
            y.append(hist.GetBinContent(i))
            xerr_down.append(hist.GetBinWidth(i) / 2.)
            xerr_up.append(hist.GetBinWidth(i) / 2.)
            yerr_down.append(hist.GetBinErrorLow(i))
            yerr_up.append(hist.GetBinErrorUp(i))

        graph = ROOT.TGraphAsymmErrors(len(x), array("f", x), array("f", y), array("f", xerr_down),
            array("f", xerr_up), array("f", yerr_down), array("f", yerr_up))
        graph.SetFillStyle(1001)
        graph.SetFillColor(16)
        graph.SetLineColor(0)
        graph.SetMarkerColor(0)

        self.add_object(graph)
        graph.Draw(" ".join(options))

    def save(self, draw_legend=False, log_y=False):
        self.pad.cd()
        if log_y:
            self.pad.SetLogy()
        if draw_legend:
            self.legend.Draw()
        for obj in self.objects:
            if hasattr(obj, "Update"):
                obj.Update()
        self.pad.Update()

    def close(self):
        self.pad.Close()

    def cd(self):
        self.pad.cd()


class ROOTPlot(object):

    def __init__(self, *args, **kwargs):
        self.canvas = ROOT.TCanvas(*args, **kwargs)

        # store all used ROOT object to make sure they are not cleaned up
        self.objects = []

        self.set_style()

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
                pad.draw_base_legend()
                self.pads[(idx_x, idx_y)] = pad
                self.open_pad = pad

    def cd(self, idx_x, idx_y):
        self.canvas.cd()
        self.open_pad = self.pads[(idx_x, idx_y)]
        self.open_pad.cd()

    def draw(self, *args, **kwargs):
        self.open_pad.draw(*args, **kwargs)

    def draw_as_graph(self, *args, **kwargs):
        self.open_pad.draw_as_graph(*args, **kwargs)

    def draw_text(self, *args, **kwargs):
        self.open_pad.draw_text(*args, **kwargs)

    def save(self, path, lumi=-1., add_cms_label=True, **kwargs):
        for pad in self.pads.values():
            pad.save(**kwargs)

        left_margin = self.canvas.GetLeftMargin()
        top_margin = self.canvas.GetTopMargin()
        right_margin = self.canvas.GetRightMargin()

        # see https://twiki.cern.ch/twiki/bin/view/CMS/Internal/FigGuidelines
        if lumi > 0:
            self.canvas.cd()
            self.lumi = ROOT.TLatex()
            self.lumi.SetTextSize(0.6*top_margin)
            self.lumi.SetTextFont(42)
            self.lumi.SetTextAlign(31)
            self.lumi.DrawLatex(1 - right_margin, 1 - 0.8*top_margin, "%.2f fb^{-1}(13 TeV)" % lumi)

        if add_cms_label:
            cms_text_size = 0.75 * top_margin
            self.cms_label = ROOT.TLatex()
            self.cms_label.SetTextSize(cms_text_size)
            self.cms_label.SetTextFont(61)
            self.cms_label.SetTextAlign(13)
            self.cms_label.DrawLatex(1.6 * left_margin, 1 - 1.2 * top_margin, "CMS")

            self.preliminary_label = ROOT.TLatex()
            self.preliminary_label.SetTextSize( 0.76 * cms_text_size)
            self.preliminary_label.SetTextFont(52)
            self.preliminary_label.SetTextAlign(13)
            self.preliminary_label.DrawLatex(1.6 * left_margin, 1 - 1.2 * top_margin - 1.5 * cms_text_size, "Preliminary")

        self.canvas.Update()
        self.canvas.SaveAs(path)

    def __del__(self):
        for pad in self.pads.values():
            pad.close()

    def set_style(self):
        # set plot style
        # adapted from https://twiki.cern.ch/twiki/bin/view/CMS/Internal/FigGuidelines
        style = ROOT.TStyle("plot style","Style for jtsf plots.")

        style.SetOptStat(0)
        style.SetOptTitle(0)
        style.SetLegendBorderSize(0)

        #for the canvas:
        style.SetCanvasBorderMode(0)
        style.SetCanvasColor(ROOT.kWhite)
        style.SetCanvasDefH(600) # Height of canvas
        style.SetCanvasDefW(600) # Width of canvas
        style.SetCanvasDefX(0) # Position on screen
        style.SetCanvasDefY(0)

        style.SetPadBorderMode(0)
        style.SetPadColor(ROOT.kWhite)
        style.SetPadGridX(False)
        style.SetPadGridY(False)
        style.SetGridColor(0)
        style.SetGridStyle(3)
        style.SetGridWidth(1)

        #For the frame:
        style.SetFrameBorderMode(0)
        style.SetFrameBorderSize(1)
        style.SetFrameFillColor(0)
        style.SetFrameFillStyle(0)
        style.SetFrameLineColor(1)
        style.SetFrameLineStyle(1)
        style.SetFrameLineWidth(1)

        #For the histo:
        style.SetHistLineColor(1)
        style.SetHistLineStyle(0)
        style.SetHistLineWidth(1)

        #style.SetEndErrorSize(2)
        #style.SetMarkerStyle(20)

        #For the fit/function:
        style.SetOptFit(1)
        style.SetFitFormat("5.4g")
        style.SetFuncColor(2)
        style.SetFuncStyle(1)
        style.SetFuncWidth(1)

        #For the date:
        style.SetOptDate(0)

        # Margins:
        style.SetPadTopMargin(0.05)
        style.SetPadBottomMargin(0.13)
        style.SetPadLeftMargin(0.16)
        style.SetPadRightMargin(0.02)

        # For the axis titles:

        style.SetTitleColor(1, "XYZ")
        style.SetTitleFont(42, "XYZ")
        style.SetTitleSize(0.06, "XYZ")
        style.SetTitleXOffset(0.9)
        style.SetTitleYOffset(1.25)

        # For the axis labels:
        style.SetLabelColor(1, "XYZ")
        style.SetLabelFont(42, "XYZ")
        style.SetLabelOffset(0.007, "XYZ")
        style.SetLabelSize(0.05, "XYZ")

        # For the axis:
        style.SetAxisColor(1, "XYZ")
        style.SetStripDecimals(True)
        style.SetTickLength(0.03, "XYZ")
        style.SetNdivisions(510, "XYZ")
        #style.SetPadTickX(1)
        #style.SetPadTickY(1)

        # Change for log plots:
        style.SetOptLogx(0)
        style.SetOptLogy(0)
        style.SetOptLogz(0)

        # Postscript options:
        style.SetPaperSize(20.,20.)

        style.SetHatchesLineWidth(5)
        style.SetHatchesSpacing(0.05)

        self.style = style
        self.style.cd()
