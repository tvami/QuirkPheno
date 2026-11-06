#!/usr/bin/env python3
"""Same-model comparison with FASER: SU(2)-singlet quirk, N_IC = 2, excluded mass vs Lambda (this recast vs FASER)."""
import json
import math
import os
import sys
import numpy as np
import ROOT
import yaml
import pstyle as P
import excl as X

NIC = 2
FARINA = "farina" in sys.argv[1:]  # add the Farina-Low singlet points (arXiv:1703.00912, Table 2)
OUT = os.environ.get("QP_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))  # outputs
res = json.load(open(OUT + "/cards/nom/limits_singlet.json"))
FASER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hepdata", "faser", "exclusion_limit_contours.yaml")


def crossing(ms, rs, thr):
    """Largest mass with r < thr (log-linear interpolation), 0 if not even the lowest mass."""
    o = np.argsort(ms)
    ms, rs = np.array(ms, float)[o], np.array(rs, float)[o]
    if rs[0] >= thr:
        return 0.
    for i in range(len(ms) - 1):
        if rs[i] < thr <= rs[i + 1]:
            return ms[i] + (ms[i + 1] - ms[i]) * (math.log(thr / rs[i])) / math.log(rs[i + 1] / rs[i])
    return ms[-1]


tab = {}
for L in sorted({r["L"] for r in res}):
    rr = [r for r in res if r["L"] == L and r["best"]]
    allm = [r["m"] for r in res if r["L"] == L]
    iv = {k: X.interval([r["m"] for r in rr], [r["best"][k] for r in rr], NIC) if len(rr) >= 2 else (0., 0.) for k in ("obs", "exp", "m1", "p1")}
    tab[L] = {k: v[1] for k, v in iv.items()}  # upper edge
    tab[L]["lo"] = iv["obs"][0]
    tab[L]["mmin"] = min(allm)
json.dump({str(k): v for k, v in tab.items()}, open(os.path.join(P.OUT, "mass_limits_singlet_nic2.json"), "w"), indent=1)
for L, t in sorted(tab.items()):
    print("Lambda %7g eV: singlet N_IC=2 obs %5.0f-" % (L, t["lo"]), "%5.0f exp %5.0f [%5.0f, %5.0f]  (lowest mass %d)" % (t["obs"], t["exp"], t["p1"], t["m1"], t["mmin"]))

Ls = sorted(tab)
MMIN = {L: tab[L]["mmin"] for L in Ls}


def segments(vals):
    """Runs with an exclusion above the lowest simulated mass, closed at that mass."""
    segs, cur = [], []
    for i, (L, v) in enumerate(zip(Ls, vals)):
        if v > MMIN[L]:
            if not cur and i > 0:
                cur.append((Ls[i - 1], MMIN[Ls[i - 1]]))
            cur.append((L, v))
        elif cur:
            cur.append((L, MMIN[L]))
            segs.append(cur)
            cur = []
    if cur:
        segs.append(cur)
    return segs


def band_poly(seg, lo, hi):
    xs = [lo[L] for L, _ in seg] + [hi[L] for L, _ in reversed(seg)]
    ys = [L for L, _ in seg] + [L for L, _ in reversed(seg)]
    return P.graph(xs, ys)


P.label("13 TeV, 101 fb^{-1}")
c = P.canvas("c_faser", 0., 1300., 0.5, 2e5, "m_{Q} [GeV]", "#Lambda [eV]", logy=True, iPos=0)
m0 = min(MMIN.values())
box = ROOT.TBox(0., 0.5, m0, 2e5)  # below the lowest mass simulated in this recast
box.SetFillColorAlpha(ROOT.kGray, 0.6)
box.SetFillStyle(1001)
box.Draw()
P.keep.append(box)
f = yaml.safe_load(open(FASER))
fl = [float(v["value"]) for v in f["independent_variables"][0]["values"]]
dv = {d["header"]["name"]: [float(v["value"]) for v in d["values"]] for d in f["dependent_variables"]}
fexp = dv["Expected excluded quirk mass at 95% CL"]
fp1 = dv["Expected excluded quirk mass at 95% CL, +1 sigma"]
fm1 = dv["Expected excluded quirk mass at 95% CL, -1 sigma"]
fobs = dv["Observed excluded quirk mass at 95% CL"]
fx = P.graph([0.] + fobs + [0.], [fl[0]] + fl + [fl[-1]])  # excluded: below the observed mass
fx.SetFillColorAlpha(ROOT.kOrange + 7, 0.25)
fx.SetFillStyle(1001)
fx.Draw("F same")
fb = P.graph(fm1 + fp1[::-1], fl + fl[::-1])
fb.SetFillColor(ROOT.kOrange + 7)
fb.SetFillStyle(3345)
fb.SetLineColor(ROOT.kOrange + 7)
fb.Draw("F same")
gfe = P.graph(fexp, fl)
gfe.SetLineColor(ROOT.kOrange + 7)
gfe.SetLineStyle(2)
gfe.SetLineWidth(2)
gfe.Draw("L same")
gfo = P.graph(fobs, fl)
gfo.SetLineColor(ROOT.kOrange + 7)
gfo.SetLineWidth(3)
gfo.Draw("L same")
import plane as PL
obs = PL.draw(res, "obs", NIC, "hs_obs", fill=ROOT.kAzure + 2, color=ROOT.kAzure + 2, width=3)
exp = PL.draw(res, "exp", NIC, "hs_exp", color=ROOT.kAzure + 2, style=2)
band = PL.draw_band(res, NIC, ROOT.kAzure + 2, "hs_band")
rx = obs[0].Clone()
rx.SetFillColorAlpha(ROOT.kAzure + 2, 0.25)
P.keep.append(rx)
go = obs[0]
gfl = None
if FARINA:
    gfl = P.graph([650., 150.], [1., 100.])  # Nc = 2 singlet, 12.9 fb^-1 of 2016 data
    gfl.SetMarkerStyle(33)
    gfl.SetMarkerSize(2.2)
    gfl.SetMarkerColor(ROOT.kGreen + 2)
    gfl.Draw("P same")
lg = P.legend(0.46, 0.54 if FARINA else 0.57, 0.93, 0.88, 0.023)
lg.AddEntry(rx, "Recast, excluded", "f")
lg.AddEntry(go, "Recast, observed", "l")
lg.AddEntry(exp[0], "Recast, expected", "l")
lg.AddEntry(band, "Recast, expected 68% band", "f")
lg.AddEntry(fx, "FASER (13.6 TeV), excluded", "f")
lg.AddEntry(gfo, "FASER (13.6 TeV), observed", "l")
lg.AddEntry(fb, "FASER (13.6 TeV), expected 68% band", "f")
if gfl:
    lg.AddEntry(gfl, "Farina, Low (2016 data), mass limit", "p")
lg.AddEntry(box, "Not simulated (recast)", "f")
lg.Draw()
P.note("SU(2)-singlet quirk, N_{IC} = %d" % NIC, 0.46, 0.49 if FARINA else 0.52, 0.028)
P.save(c, "faser_farina_comparison" if FARINA else "faser_comparison")
