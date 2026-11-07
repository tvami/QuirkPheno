#!/usr/bin/env python3
"""Quirk limits: r = sigma95/sigma_LO vs mass per Lambda, and the excluded mass vs Lambda for N_IC = 1, 2, 3
(doublet; the FASER comparison is in faser_cmp.py). usage: limits.py [cards tag]"""
import json
import math
import os
import sys
import numpy as np
import ROOT
import pstyle
import excl

TAG = sys.argv[1] if len(sys.argv) > 1 else "nom"
OUT = os.environ.get("QP_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))  # outputs
res = json.load(open(OUT + "/cards/%s/limits.json" % TAG))


def crossing(ms, rs, thr):
    """Largest mass with r < thr, log-linear interpolation; 0 if never, inf-like if always."""
    ms, rs = np.array(ms, float), np.array(rs, float)
    order = np.argsort(ms)
    ms, rs = ms[order], rs[order]
    if rs[0] >= thr:
        return 0.
    for idx in range(len(ms) - 1):
        if rs[idx] < thr <= rs[idx + 1]:
            log0, log1 = math.log(rs[idx]), math.log(rs[idx + 1])
            return ms[idx] + (ms[idx + 1] - ms[idx]) * (math.log(thr) - log0) / (log1 - log0)
    return ms[-1]  # excluded up to the last grid point


lams = sorted({row["L"] for row in res})
tab = {}
for lam in lams:
    rr = sorted([row for row in res if row["L"] == lam and row["best"]], key=lambda row: row["m"])
    if len(rr) < 2:  # no signal passes anywhere: nothing excluded
        tab[lam] = {key: {nic: 0. for nic in (1, 2, 3)} for key in ("obs", "exp", "m1", "p1")}
        tab[lam]["lo"] = {key: {nic: 0. for nic in (1, 2, 3)} for key in ("obs", "exp", "m1", "p1")}
        tab[lam].update(mmax=max(row["m"] for row in res if row["L"] == lam), methods=[])
        continue
    ms = [row["m"] for row in rr]
    iv = {key: {nic: excl.interval(ms, [row["best"][key] for row in rr], nic) for nic in (1, 2, 3)}
          for key in ("obs", "exp", "m1", "p1")}
    tab[lam] = {key: {nic: iv[key][nic][1] for nic in (1, 2, 3)} for key in iv}  # upper edge
    tab[lam]["lo"] = {key: {nic: iv[key][nic][0] for nic in (1, 2, 3)} for key in iv}  # lower edge
    tab[lam]["mmax"] = max(ms)
    tab[lam]["methods"] = [row["method"] for row in rr]
json.dump({str(key): val for key, val in tab.items()}, open(os.path.join(pstyle.OUT, "mass_limits_%s.json" % TAG), "w"), indent=1)
for lam, lim in sorted(tab.items()):
    print("Lambda %7g eV: N_IC=1 obs %5.0f-%5.0f exp %5.0f [%5.0f, %5.0f]  N_IC=2 obs %5.0f  N_IC=3 obs %5.0f  (grid to %d)" % (
        lam, lim["lo"]["obs"][1], lim["obs"][1], lim["exp"][1], lim["p1"][1], lim["m1"][1],
        lim["obs"][2], lim["obs"][3], lim["mmax"]))

# r vs mass per Lambda
pstyle.label("13 TeV, 101 fb^{-1}")
canv = pstyle.canvas("c_r", 0, 1500, 1e-3, 1e6, "m_{Q} [GeV]", "#sigma_{95} / #sigma_{LO} (N_{IC} = 1)", logy=True, iPos=0)
lg = pstyle.legend(0.24, 0.70, 0.93, 0.88, 0.028)
lg.SetNColumns(2)
SHOW = [1., 30., 50., 70., 100., 10000., 30000.]
for idx, lam in enumerate([lam for lam in lams if lam in SHOW]):
    rr = sorted([row for row in res if row["L"] == lam and row["best"]], key=lambda row: row["m"])
    if not rr:
        continue
    gr = pstyle.graph([row["m"] for row in rr], [min(max(row["best"]["obs"], 1.1e-3), 9e5) for row in rr])
    col = [ROOT.kBlack, ROOT.kAzure + 1, ROOT.kOrange + 7, ROOT.kGreen + 2, ROOT.kRed + 1, ROOT.kViolet + 1, ROOT.kCyan + 2][idx % 7]
    gr.SetLineColor(col)
    gr.SetMarkerColor(col)
    gr.SetMarkerStyle([20, 21, 22, 23, 33, 29, 34][idx % 7])
    gr.SetLineWidth(2)
    gr.Draw("LP same")
    lg.AddEntry(gr, "#Lambda = %s" % ("%g eV" % lam if lam < 1000 else "%g keV" % (lam / 1000.)), "lp")
ln = ROOT.TLine(0, 1., 1500, 1.)
ln.SetLineStyle(2)
ln.Draw()
pstyle.keep.append(ln)
lg.Draw()
pstyle.save(canv, "limit_r_vs_mass")

# excluded region in the (m_Q, Lambda) plane: smooth contours of r = N_IC
import plane
MMIN = min(row["m"] for row in res)
c2 = pstyle.canvas("c_ml", 0., 1500., 0.5, 2e5, "m_{Q} [GeV]", "#Lambda [eV]", logy=True, iPos=0)
box = ROOT.TBox(0., 0.5, MMIN, 2e5)
box.SetFillColorAlpha(ROOT.kGray, 0.6)
box.SetFillStyle(1001)
box.Draw()
pstyle.keep.append(box)
obs = plane.draw(res, "obs", 1., "h_obs1", fill=ROOT.kAzure + 2, width=3)
exp = plane.draw(res, "exp", 1., "h_exp1", color=ROOT.kGreen + 2, style=2)
band = plane.draw_band(res, 1., ROOT.kGreen + 2, "h_band")
o2 = plane.draw(res, "obs", 2., "h_obs2", style=7)
o3 = plane.draw(res, "obs", 3., "h_obs3", style=3)
fx = obs[0].Clone()
fx.SetFillColorAlpha(ROOT.kAzure + 2, 0.25)
pstyle.keep.append(fx)
lg2 = pstyle.legend(0.50, 0.62, 0.93, 0.88, 0.028)
lg2.AddEntry(fx, "Excluded, N_{IC} = 1", "f")
lg2.AddEntry(obs[0], "Observed, N_{IC} = 1", "l")
lg2.AddEntry(exp[0], "Expected, N_{IC} = 1", "l")
lg2.AddEntry(band, "Expected 68% band", "f")
lg2.AddEntry(o2[0] if o2 else obs[0], "Observed, N_{IC} = 2", "l")
lg2.AddEntry(o3[0] if o3 else obs[0], "Observed, N_{IC} = 3", "l")
lg2.Draw()
pstyle.save(c2, "limit_mass_vs_lambda")
