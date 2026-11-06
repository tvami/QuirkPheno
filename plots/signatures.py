#!/usr/bin/env python3
"""Quirk signatures vs Lambda: trigger, preselection and signal-region efficiencies, merged-track fraction, Gstrip."""
import glob
import os
import re
import sys
import numpy as np
import ROOT
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "recast"))
import pstyle
import summarize

OUT = os.environ.get("QP_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))  # outputs
EMU = OUT + "/emu/"
pts = {}
for fname in (sorted(glob.glob(EMU + "quirk_m*_L*.npz"))
              + sorted(glob.glob(EMU + "quirkb_m*_L*.npz"))):  # quirkb (bound filter) overrides
    mass, lam = map(float, re.search(r"_m(\d+)_L([\d.]+)\.npz", fname).groups())
    npz = np.load(fname)
    if int(npz["nsim"]) == 0:
        continue
    col = {name: pos for pos, name in enumerate(npz["keys"])}
    arr = npz["data"]
    ng = float(npz["ngen"])
    trig = arr[:, col["trig"]]
    pre = trig * (arr[:, col["depth"]] >= summarize.NSTEP)
    ion = sum(vals.sum() for vals in summarize.bins(fname).values())
    cm = (arr[:, col["depth"]] >= summarize.NSTEP)
    pts[(int(mass), lam)] = dict(trig=trig.sum() / ng, pre=pre.sum() / ng, ion=ion,
                                 merged=float((arr[cm, col["c_merged"]]).mean()) if cm.any() else 0., f=fname)
lams = sorted({pt[1] for pt in pts})
print("points:", len(pts), "lambdas:", lams)
pstyle.label("13 TeV")
for key, yt, stem in (("trig", "Trigger efficiency", "sig_trig"), ("pre", "Trigger and preselection efficiency", "sig_presel"),
                      ("ion", "Signal efficiency, ionization method", "sig_ion")):
    canv = pstyle.canvas("c_" + key, 0.5, 5e4, 1e-4, 2., "#Lambda [eV]", yt, logy=True, logx=True, iPos=0)
    lg = pstyle.legend(0.65, 0.70, 0.93, 0.88, 0.032)
    for idx, mass in enumerate((200, 600, 1000, 1400)):
        xs = [lam for lam in lams if (mass, lam) in pts]
        if not xs:
            continue
        gr = pstyle.graph(xs, [max(pts[(mass, lam)][key], 1.1e-4) for lam in xs])
        gr.SetLineColor(pstyle.COLORS[idx])
        gr.SetMarkerColor(pstyle.COLORS[idx])
        gr.SetLineWidth(2)
        gr.SetMarkerStyle(20)
        gr.Draw("LP same")
        lg.AddEntry(gr, "m_{Q} = %d GeV" % mass, "lp")
    lg.Draw()
    pstyle.save(canv, stem)
for mass in (200, 600, 1000, 1400):
    print(mass, " ".join("L%g:trig %.3f pre %.3f ion %.3f merged %.2f" % (
        lam, pts[(mass, lam)]["trig"], pts[(mass, lam)]["pre"], pts[(mass, lam)]["ion"], pts[(mass, lam)]["merged"])
        for lam in lams if (mass, lam) in pts))

# Gstrip of the selected candidates at 1 TeV
canv = pstyle.canvas("c_gs", 0., 1., 0., 0.5, "G_{i}^{Strips}", "Fraction of selected candidates", iPos=0)
lg = pstyle.legend(0.20, 0.66, 0.60, 0.88, 0.032)
for idx, lam in enumerate([lam for lam in (1., 100., 1000., 10000.) if (1000, lam) in pts]):
    npz = np.load(pts[(1000, lam)]["f"])
    col = {name: pos for pos, name in enumerate(npz["keys"])}
    arr = npz["data"]
    sel = (arr[:, col["depth"]] >= summarize.NSTEP) & (arr[:, col["c_pt"]] > 200)
    if sel.sum() == 0:
        continue
    histo = pstyle.hist("hgs%d" % idx, arr[sel, col["c_gstrip"]], 20, 0., 1., arr[sel, col["trig"]])
    histo.Scale(1. / max(histo.Integral(), 1e-9))
    histo.SetLineColor(pstyle.COLORS[idx])
    histo.SetLineWidth(2)
    histo.Draw("hist same")
    lg.AddEntry(histo, "#Lambda = %g eV" % lam, "l")
lg.Draw()
pstyle.note("m_{Q} = 1 TeV, p_{T} > 200 GeV", 0.20, 0.60)
pstyle.save(canv, "sig_gstrip")
