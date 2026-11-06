#!/usr/bin/env python3
"""Transverse view of the gun quirk pair (m 200 GeV) through the CMS-like detector, all quirk steps."""
import math
import numpy as np
import ROOT
import uproot
import pstyle as P

D = "/ceph/cms/store/user/tvami/HSCP-Quirks/QuirkPheno/display/"
RADII = [1.16, 1.29, 1.77, 2.95, 4.2, 7.4]  # tracker, ECAL, HCAL, solenoid, muon stations
P.label("CMS-like detector, Geant4")
for L in (30, 100, 300):
    s = uproot.open(D + "ed_L%d.root" % L)["steps"].arrays(["q", "x", "y", "t"], library="np")
    c = P.canvas("c_ed%d" % L, -1., 8., -4.5, 4.5, "x [m]", "y [m]", iPos=0)
    for r in RADII:
        pts = [(r * math.cos(2 * math.pi * k / 720), r * math.sin(2 * math.pi * k / 720)) for k in range(721)]
        arcs, cur = [], []
        for p in pts:  # split the circle where it leaves the frame
            if -1 < p[0] < 8 and -4.5 < p[1] < 4.5:
                cur.append(p)
            elif cur:
                arcs.append(cur)
                cur = []
        arcs.append(cur)
        for a in arcs:
            if len(a) < 2:
                continue
            e = P.graph([p[0] for p in a], [p[1] for p in a])
            e.SetLineColor(ROOT.kGray + 1)
            e.SetLineStyle(2)
            e.Draw("L same")
    lg = P.legend(0.62, 0.20, 0.93, 0.30, 0.03)
    for q, col, name in ((0, ROOT.kAzure + 1, "quirk"), (1, ROOT.kOrange + 7, "antiquirk")):
        m = s["q"] == q
        o = np.argsort(s["t"][m])
        g = P.graph(s["x"][m][o] / 1000., s["y"][m][o] / 1000.)
        g.SetLineColor(col)
        g.SetLineWidth(2)
        g.Draw("L same")
        lg.AddEntry(g, name, "l")
    lg.Draw()
    P.note("m_{Q} = 200 GeV, #Lambda = %d eV" % L, 0.20, 0.85)
    P.save(c, "ed_xy_L%d" % L)
