#!/usr/bin/env python3
"""MIP dE/dx templates per sensor thickness from the 10 GeV pion-mass gun (mip_eta*.root), noise added."""
import glob
import os
import numpy as np
import qreco as Q

OUT = os.environ.get("QP_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))  # outputs
D = os.path.join(OUT, "mip")
rng = np.random.default_rng(11)
meta, thick = Q.load_meta(D + "/mip_eta0.05_meta.json")
vals = {}
for f in sorted(glob.glob(D + "/mip_eta*.root")):
    for a in Q.read_events(f):
        for i in range(len(a["evt"])):
            det, vol, dE, ln = a["h_det"][i], a["h_vol"][i], a["h_dE"][i], a["h_len"][i]
            m = np.isin(det, Q.TRACKER) & (ln > 0)
            for v, e, l, dd in zip(vol[m], dE[m], ln[m], det[m]):
                x = e / (l / 10.)  # MeV/cm
                x *= 1 + rng.normal(0., Q.NOISE["pix" if dd in Q.PIXEL else "strip"])
                vals.setdefault(int(thick[v]), []).append(x)
out = {"t%d" % k: np.sort(np.array(v)) for k, v in vals.items()}
np.savez(Q.HERE + "/mip_templates_%s.npz" % Q.TAG, **out)
for k, v in sorted(out.items()):
    print(k, len(v), "median %.3f MeV/cm  q10 %.3f q90 %.3f" % (np.median(v), np.quantile(v, .1), np.quantile(v, .9)))
