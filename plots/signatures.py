#!/usr/bin/env python3
"""Quirk signatures vs Lambda: trigger, preselection and signal-region efficiencies, merged-track fraction, Gstrip."""
import glob
import os
import re
import sys
import numpy as np
import ROOT
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "recast"))
import pstyle as P
import summarize as S

E = "/ceph/cms/store/user/tvami/HSCP-Quirks/QuirkPheno/emu/"
pts = {}
for f in sorted(glob.glob(E + "quirk_m*_L*.npz")) + sorted(glob.glob(E + "quirkb_m*_L*.npz")):  # quirkb (bound filter) overrides
    m, L = map(float, re.search(r"_m(\d+)_L([\d.]+)\.npz", f).groups())
    d = np.load(f)
    if int(d["nsim"]) == 0:
        continue
    k = {s: i for i, s in enumerate(d["keys"])}
    x = d["data"]
    ng = float(d["ngen"])
    trig = x[:, k["trig"]]
    pre = trig * (x[:, k["depth"]] >= S.NSTEP)
    ion = sum(v.sum() for v in S.bins(f).values())
    cm = (x[:, k["depth"]] >= S.NSTEP)
    pts[(int(m), L)] = dict(trig=trig.sum() / ng, pre=pre.sum() / ng, ion=ion,
                            merged=float((x[cm, k["c_merged"]]).mean()) if cm.any() else 0., f=f)
lams = sorted({k[1] for k in pts})
print("points:", len(pts), "lambdas:", lams)
P.label("13 TeV")
for key, yt, stem in (("trig", "Trigger efficiency", "sig_trig"), ("pre", "Trigger and preselection efficiency", "sig_presel"),
                      ("ion", "Signal efficiency, ionization method", "sig_ion")):
    c = P.canvas("c_" + key, 0.5, 5e4, 1e-4, 2., "#Lambda [eV]", yt, logy=True, logx=True, iPos=0)
    lg = P.legend(0.65, 0.70, 0.93, 0.88, 0.032)
    for i, m in enumerate((200, 600, 1000, 1400)):
        xs = [L for L in lams if (m, L) in pts]
        if not xs:
            continue
        g = P.graph(xs, [max(pts[(m, L)][key], 1.1e-4) for L in xs])
        g.SetLineColor(P.COLORS[i])
        g.SetMarkerColor(P.COLORS[i])
        g.SetLineWidth(2)
        g.SetMarkerStyle(20)
        g.Draw("LP same")
        lg.AddEntry(g, "m_{Q} = %d GeV" % m, "lp")
    lg.Draw()
    P.save(c, stem)
for m in (200, 600, 1000, 1400):
    print(m, " ".join("L%g:trig %.3f pre %.3f ion %.3f merged %.2f" % (L, pts[(m, L)]["trig"], pts[(m, L)]["pre"], pts[(m, L)]["ion"], pts[(m, L)]["merged"]) for L in lams if (m, L) in pts))

# Gstrip of the selected candidates at 1 TeV
c = P.canvas("c_gs", 0., 1., 0., 0.5, "G_{i}^{Strips}", "Fraction of selected candidates", iPos=0)
lg = P.legend(0.20, 0.66, 0.60, 0.88, 0.032)
for i, L in enumerate([L for L in (1., 100., 1000., 10000.) if (1000, L) in pts]):
    d = np.load(pts[(1000, L)]["f"])
    k = {s: j for j, s in enumerate(d["keys"])}
    x = d["data"]
    s = (x[:, k["depth"]] >= S.NSTEP) & (x[:, k["c_pt"]] > 200)
    if s.sum() == 0:
        continue
    h = P.hist("hgs%d" % i, x[s, k["c_gstrip"]], 20, 0., 1., x[s, k["trig"]])
    h.Scale(1. / max(h.Integral(), 1e-9))
    h.SetLineColor(P.COLORS[i])
    h.SetLineWidth(2)
    h.Draw("hist same")
    lg.AddEntry(h, "#Lambda = %g eV" % L, "l")
lg.Draw()
P.note("m_{Q} = 1 TeV, p_{T} > 200 GeV", 0.20, 0.60)
P.save(c, "sig_gstrip")
