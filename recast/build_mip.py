#!/usr/bin/env python3
"""MIP dE/dx templates per sensor thickness from the 10 GeV pion-mass gun (mip_eta*.root), noise added."""
import glob
import os
import numpy as np
import qreco

OUT = os.environ.get("QP_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))  # outputs
MIPDIR = os.path.join(OUT, "mip")
rng = np.random.default_rng(11)
meta, thick = qreco.load_meta(MIPDIR + "/mip_eta0.05_meta.json")
vals = {}
for fname in sorted(glob.glob(MIPDIR + "/mip_eta*.root")):
    for arrays in qreco.read_events(fname):
        for ievt in range(len(arrays["evt"])):
            det, vol, dE, ln = arrays["h_det"][ievt], arrays["h_vol"][ievt], arrays["h_dE"][ievt], arrays["h_len"][ievt]
            mask = np.isin(det, qreco.TRACKER) & (ln > 0)
            for vid, edep, length, dd in zip(vol[mask], dE[mask], ln[mask], det[mask]):
                dedx = edep / (length / 10.)  # MeV/cm
                dedx *= 1 + rng.normal(0., qreco.NOISE["pix" if dd in qreco.PIXEL else "strip"])
                vals.setdefault(int(thick[vid]), []).append(dedx)
out = {"t%d" % key: np.sort(np.array(tmpl)) for key, tmpl in vals.items()}
np.savez(qreco.HERE + "/mip_templates_%s.npz" % qreco.TAG, **out)
for key, tmpl in sorted(out.items()):
    print(key, len(tmpl), "median %.3f MeV/cm  q10 %.3f q90 %.3f"
          % (np.median(tmpl), np.quantile(tmpl, .1), np.quantile(tmpl, .9)))
