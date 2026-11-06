#!/usr/bin/env python3
"""Gstrip in the PASS region (pT > 200 GeV): CMS data and predicted background (HEPData Fig. 4 right) with emulated quirk signals."""
import array
import os
import sys
import numpy as np
import ROOT
import yaml
import pstyle as P
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "recast"))
import summarize as S

B = "/ceph/cms/store/user/tvami/HSCP-Quirks/QuirkPheno/"
LUMI = 101.0e3  # pb^-1
d = yaml.safe_load(open(os.path.join(S.Q.HEPDATA, "figure_4_right.yaml")))
edges = [float(v["low"]) for v in d["independent_variables"][0]["values"]] + [float(d["independent_variables"][0]["values"][-1]["high"])]
dv = {x["header"]["name"]: [float(v["value"]) for v in x["values"]] for x in d["dependent_variables"]}
data, bkg = dv["Number of data events"], dv["Number of background events"]
SIGS = [("quirk", 800, 1, "m_{Q} = 800 GeV, #Lambda = 1 eV", ROOT.kAzure + 1, 1),
        ("quirk", 600, 50, "m_{Q} = 600 GeV, #Lambda = 50 eV", ROOT.kGreen + 2, 7),
        ("quirkb", 200, 10000, "m_{Q} = 200 GeV, #Lambda = 10 keV", ROOT.kOrange + 7, 1),
        ("quirkb", 300, 30000, "m_{Q} = 300 GeV, #Lambda = 30 keV", ROOT.kViolet + 1, 2)]
ed = array.array("d", edges)


def hist(name, vals):
    h = ROOT.TH1D(name, "", len(vals), ed)
    h.SetDirectory(0)
    for i, v in enumerate(vals):
        h.SetBinContent(i + 1, v)
    P.keep.append(h)
    return h


P.label("13 TeV, 101 fb^{-1}")
c = P.canvas("c_gs_pass", 0., 1., 0.05, 5e6, "G_{i}^{Strips}", "Events / bin", logy=True, iPos=0)
hb = hist("hbkg", bkg)
hb.SetFillColor(ROOT.kGray)
hb.SetLineColor(ROOT.kGray + 2)
hb.Draw("hist same")
lg = P.legend(0.40, 0.62, 0.93, 0.88, 0.028)
lg.AddEntry(hb, "Predicted background, CMS [6]", "f")
for model, m, L, lab, col, sty in SIGS:
    fn = B + "emu/%s_m%d_L%d.npz" % (model, m, L)
    xs = B + "cards/nom/xs_%s_m%d_L%d.txt" % (model, m, L)
    if not (os.path.exists(fn) and os.path.exists(xs)):
        print("missing", fn)
        continue
    sigma = float(open(xs).read().split()[2])
    b = S.bins(fn)
    hs = hist("hs_%d_%d" % (m, L), list(np.concatenate([b["pass_LOW"], b["pass_SIG"]]) * sigma * LUMI))
    hs.SetLineColor(col)
    hs.SetLineStyle(sty)
    hs.SetLineWidth(3)
    hs.Draw("hist same")
    lg.AddEntry(hs, lab, "l")
    print("%-34s signal per bin" % lab, np.round(np.concatenate([b["pass_LOW"], b["pass_SIG"]]) * sigma * LUMI, 2))
gd = ROOT.TGraphAsymmErrors(len(data))
for i, n in enumerate(data):
    gd.SetPoint(i, 0.5 * (edges[i] + edges[i + 1]), n)
    gd.SetPointError(i, 0, 0, ROOT.Math.gamma_quantile(0.1587, n, 1.) if n > 0 else 0., 0)
    lo = n - ROOT.Math.gamma_quantile(0.1587, n, 1.) if n > 0 else 0.
    hi = ROOT.Math.gamma_quantile_c(0.1587, n + 1, 1.) - n  # Garwood intervals
    gd.SetPointError(i, 0, 0, lo, hi)
gd.SetMarkerStyle(20)
gd.Draw("P same")
P.keep.append(gd)
lg.AddEntry(gd, "Data, CMS [6]", "pe")  # ref number of CMS:EXO18002 in the paper
lg.Draw()
P.note("PASS (F_{i}^{Pixels} > 0.9), p_{T} > 200 GeV", 0.40, 0.57, 0.028)
P.note("Quirks: SU(2) doublet, N_{IC} = 1, LO", 0.40, 0.53, 0.028)
P.save(c, "gstrip_pass")
