#!/usr/bin/env python3
"""Gstrip in the PASS region (pT > 200 GeV): CMS data and predicted background (HEPData Fig. 4 right) with emulated quirk signals."""
import array
import os
import sys
import numpy as np
import ROOT
import yaml
import pstyle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "recast"))
import summarize
import qreco

OUT = os.environ.get("QP_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))  # outputs
BASE = OUT + "/"
LUMI = 101.0e3  # pb^-1
fig = yaml.safe_load(open(os.path.join(qreco.HEPDATA, "figure_4_right.yaml")))
edges = [float(val["low"]) for val in fig["independent_variables"][0]["values"]] + [float(fig["independent_variables"][0]["values"][-1]["high"])]
dv = {dep["header"]["name"]: [float(val["value"]) for val in dep["values"]] for dep in fig["dependent_variables"]}
data, bkg = dv["Number of data events"], dv["Number of background events"]
SIGS = [("quirk", 800, 1, "m_{Q} = 800 GeV, #Lambda = 1 eV", ROOT.kAzure + 1, 1),
        ("quirk", 600, 50, "m_{Q} = 600 GeV, #Lambda = 50 eV", ROOT.kGreen + 2, 7),
        ("quirkb", 200, 10000, "m_{Q} = 200 GeV, #Lambda = 10 keV", ROOT.kOrange + 7, 1),
        ("quirkb", 300, 30000, "m_{Q} = 300 GeV, #Lambda = 30 keV", ROOT.kViolet + 1, 2)]
ed = array.array("d", edges)


def hist(name, vals):
    histo = ROOT.TH1D(name, "", len(vals), ed)
    histo.SetDirectory(0)
    for idx, val in enumerate(vals):
        histo.SetBinContent(idx + 1, val)
    pstyle.keep.append(histo)
    return histo


pstyle.label("13 TeV, 101 fb^{-1}")
canv = pstyle.canvas("c_gs_pass", 0., 1., 0.05, 5e6, "G_{i}^{Strips}", "Events / bin", logy=True, iPos=0)
hb = hist("hbkg", bkg)
hb.SetFillColor(ROOT.kGray)
hb.SetLineColor(ROOT.kGray + 2)
hb.Draw("hist same")
lg = pstyle.legend(0.40, 0.62, 0.93, 0.88, 0.028)
lg.AddEntry(hb, "Predicted background, CMS [6]", "f")
for model, mass, lam, lab, col, sty in SIGS:
    fn = BASE + "emu/%s_m%d_L%d.npz" % (model, mass, lam)
    xs = BASE + "cards/nom/xs_%s_m%d_L%d.txt" % (model, mass, lam)
    if not (os.path.exists(fn) and os.path.exists(xs)):
        print("missing", fn)
        continue
    sigma = float(open(xs).read().split()[2])
    effs = summarize.bins(fn)
    hs = hist("hs_%d_%d" % (mass, lam), list(np.concatenate([effs["pass_LOW"], effs["pass_SIG"]]) * sigma * LUMI))
    hs.SetLineColor(col)
    hs.SetLineStyle(sty)
    hs.SetLineWidth(3)
    hs.Draw("hist same")
    lg.AddEntry(hs, lab, "l")
    print("%-34s signal per bin" % lab, np.round(np.concatenate([effs["pass_LOW"], effs["pass_SIG"]]) * sigma * LUMI, 2))
gd = ROOT.TGraphAsymmErrors(len(data))
for idx, nobs in enumerate(data):
    gd.SetPoint(idx, 0.5 * (edges[idx] + edges[idx + 1]), nobs)
    gd.SetPointError(idx, 0, 0, ROOT.Math.gamma_quantile(0.1587, nobs, 1.) if nobs > 0 else 0., 0)
    lo = nobs - ROOT.Math.gamma_quantile(0.1587, nobs, 1.) if nobs > 0 else 0.
    hi = ROOT.Math.gamma_quantile_c(0.1587, nobs + 1, 1.) - nobs  # Garwood intervals
    gd.SetPointError(idx, 0, 0, lo, hi)
gd.SetMarkerStyle(20)
gd.Draw("P same")
pstyle.keep.append(gd)
lg.AddEntry(gd, "Data, CMS [6]", "pe")  # ref number of CMS:EXO18002 in the paper
lg.Draw()
pstyle.note("PASS (F_{i}^{Pixels} > 0.9), p_{T} > 200 GeV", 0.40, 0.57, 0.028)
pstyle.note("Quirks: SU(2) doublet, N_{IC} = 1, LO", 0.40, 0.53, 0.028)
pstyle.save(canv, "gstrip_pass")
