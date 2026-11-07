#!/usr/bin/env python3
"""Transverse view of the gun quirk pair (m 200 GeV) through the CMS-like detector, all quirk steps."""
import math
import os
import numpy as np
import ROOT
import uproot
import pstyle

OUT = os.environ.get("QP_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))  # outputs
DISPLAY = OUT + "/display/"
RADII = [1.16, 1.29, 1.77, 2.95, 4.2, 7.4]  # tracker, ECAL, HCAL, solenoid, muon stations
pstyle.label("CMS-like detector, Geant4")
for lam in (30, 100, 300):
    steps = uproot.open(DISPLAY + "ed_L%d.root" % lam)["steps"].arrays(["q", "x", "y", "t"], library="np")
    canv = pstyle.canvas("c_ed%d" % lam, -1., 8., -4.5, 4.5, "x [m]", "y [m]", iPos=0)
    for rad in RADII:
        pts = [(rad * math.cos(2 * math.pi * ipt / 720), rad * math.sin(2 * math.pi * ipt / 720)) for ipt in range(721)]
        arcs, cur = [], []
        for pt in pts:  # split the circle where it leaves the frame
            if -1 < pt[0] < 8 and -4.5 < pt[1] < 4.5:
                cur.append(pt)
            elif cur:
                arcs.append(cur)
                cur = []
        arcs.append(cur)
        for arc in arcs:
            if len(arc) < 2:
                continue
            ring = pstyle.graph([pt[0] for pt in arc], [pt[1] for pt in arc])
            ring.SetLineColor(ROOT.kGray + 1)
            ring.SetLineStyle(2)
            ring.Draw("L same")
    lg = pstyle.legend(0.62, 0.20, 0.93, 0.30, 0.03)
    for charge, col, name in ((0, ROOT.kAzure + 1, "quirk"), (1, ROOT.kOrange + 7, "antiquirk")):
        mask = steps["q"] == charge
        order = np.argsort(steps["t"][mask])
        gr = pstyle.graph(steps["x"][mask][order] / 1000., steps["y"][mask][order] / 1000.)
        gr.SetLineColor(col)
        gr.SetLineWidth(2)
        gr.Draw("L same")
        lg.AddEntry(gr, name, "l")
    lg.Draw()
    pstyle.note("m_{Q} = 200 GeV, #Lambda = %d eV" % lam, 0.20, 0.85)
    pstyle.save(canv, "ed_xy_L%d" % lam)
