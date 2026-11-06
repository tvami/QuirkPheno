#!/usr/bin/env python3
"""Quirk limits: r = sigma95/sigma_LO vs mass per Lambda, and the excluded mass vs Lambda for N_IC = 1, 2, 3
(doublet; the FASER comparison is in faser_cmp.py). usage: limits.py [cards tag]"""
import json
import math
import os
import sys
import numpy as np
import ROOT
import pstyle as P
import excl as X

TAG = sys.argv[1] if len(sys.argv) > 1 else "nom"
OUT = os.environ.get("QP_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))  # outputs
res = json.load(open(OUT + "/cards/%s/limits.json" % TAG))


def crossing(ms, rs, thr):
    """Largest mass with r < thr, log-linear interpolation; 0 if never, inf-like if always."""
    ms, rs = np.array(ms, float), np.array(rs, float)
    o = np.argsort(ms)
    ms, rs = ms[o], rs[o]
    if rs[0] >= thr:
        return 0.
    for i in range(len(ms) - 1):
        if rs[i] < thr <= rs[i + 1]:
            a, b = math.log(rs[i]), math.log(rs[i + 1])
            return ms[i] + (ms[i + 1] - ms[i]) * (math.log(thr) - a) / (b - a)
    return ms[-1]  # excluded up to the last grid point


lams = sorted({r["L"] for r in res})
tab = {}
for L in lams:
    rr = sorted([r for r in res if r["L"] == L and r["best"]], key=lambda r: r["m"])
    if len(rr) < 2:  # no signal passes anywhere: nothing excluded
        tab[L] = {k: {n: 0. for n in (1, 2, 3)} for k in ("obs", "exp", "m1", "p1")}
        tab[L]["lo"] = {k: {n: 0. for n in (1, 2, 3)} for k in ("obs", "exp", "m1", "p1")}
        tab[L].update(mmax=max(r["m"] for r in res if r["L"] == L), methods=[])
        continue
    ms = [r["m"] for r in rr]
    iv = {k: {n: X.interval(ms, [r["best"][k] for r in rr], n) for n in (1, 2, 3)} for k in ("obs", "exp", "m1", "p1")}
    tab[L] = {k: {n: iv[k][n][1] for n in (1, 2, 3)} for k in iv}  # upper edge
    tab[L]["lo"] = {k: {n: iv[k][n][0] for n in (1, 2, 3)} for k in iv}  # lower edge
    tab[L]["mmax"] = max(ms)
    tab[L]["methods"] = [r["method"] for r in rr]
json.dump({str(k): v for k, v in tab.items()}, open(os.path.join(P.OUT, "mass_limits_%s.json" % TAG), "w"), indent=1)
for L, t in sorted(tab.items()):
    print("Lambda %7g eV: N_IC=1 obs %5.0f-%5.0f exp %5.0f [%5.0f, %5.0f]  N_IC=2 obs %5.0f  N_IC=3 obs %5.0f  (grid to %d)" % (
        L, t["lo"]["obs"][1], t["obs"][1], t["exp"][1], t["p1"][1], t["m1"][1], t["obs"][2], t["obs"][3], t["mmax"]))

# r vs mass per Lambda
P.label("13 TeV, 101 fb^{-1}")
c = P.canvas("c_r", 0, 1500, 1e-3, 1e6, "m_{Q} [GeV]", "#sigma_{95} / #sigma_{LO} (N_{IC} = 1)", logy=True, iPos=0)
lg = P.legend(0.24, 0.70, 0.93, 0.88, 0.028)
lg.SetNColumns(2)
SHOW = [1., 30., 50., 70., 100., 10000., 30000.]
for i, L in enumerate([L for L in lams if L in SHOW]):
    rr = sorted([r for r in res if r["L"] == L and r["best"]], key=lambda r: r["m"])
    if not rr:
        continue
    g = P.graph([r["m"] for r in rr], [min(max(r["best"]["obs"], 1.1e-3), 9e5) for r in rr])
    col = [ROOT.kBlack, ROOT.kAzure + 1, ROOT.kOrange + 7, ROOT.kGreen + 2, ROOT.kRed + 1, ROOT.kViolet + 1, ROOT.kCyan + 2][i % 7]
    g.SetLineColor(col)
    g.SetMarkerColor(col)
    g.SetMarkerStyle([20, 21, 22, 23, 33, 29, 34][i % 7])
    g.SetLineWidth(2)
    g.Draw("LP same")
    lg.AddEntry(g, "#Lambda = %s" % ("%g eV" % L if L < 1000 else "%g keV" % (L / 1000.)), "lp")
ln = ROOT.TLine(0, 1., 1500, 1.)
ln.SetLineStyle(2)
ln.Draw()
P.keep.append(ln)
lg.Draw()
P.save(c, "limit_r_vs_mass")

# excluded region in the (m_Q, Lambda) plane: smooth contours of r = N_IC
import plane as PL
MMIN = min(r["m"] for r in res)
c2 = P.canvas("c_ml", 0., 1500., 0.5, 2e5, "m_{Q} [GeV]", "#Lambda [eV]", logy=True, iPos=0)
box = ROOT.TBox(0., 0.5, MMIN, 2e5)
box.SetFillColorAlpha(ROOT.kGray, 0.6)
box.SetFillStyle(1001)
box.Draw()
P.keep.append(box)
obs = PL.draw(res, "obs", 1., "h_obs1", fill=ROOT.kAzure + 2, width=3)
exp = PL.draw(res, "exp", 1., "h_exp1", color=ROOT.kGreen + 2, style=2)
band = PL.draw_band(res, 1., ROOT.kGreen + 2, "h_band")
o2 = PL.draw(res, "obs", 2., "h_obs2", style=7)
o3 = PL.draw(res, "obs", 3., "h_obs3", style=3)
fx = obs[0].Clone()
fx.SetFillColorAlpha(ROOT.kAzure + 2, 0.25)
P.keep.append(fx)
lg2 = P.legend(0.50, 0.62, 0.93, 0.88, 0.028)
lg2.AddEntry(fx, "Excluded, N_{IC} = 1", "f")
lg2.AddEntry(obs[0], "Observed, N_{IC} = 1", "l")
lg2.AddEntry(exp[0], "Expected, N_{IC} = 1", "l")
lg2.AddEntry(band, "Expected 68% band", "f")
lg2.AddEntry(o2[0] if o2 else obs[0], "Observed, N_{IC} = 2", "l")
lg2.AddEntry(o3[0] if o3 else obs[0], "Observed, N_{IC} = 3", "l")
lg2.Draw()
P.save(c2, "limit_mass_vs_lambda")
