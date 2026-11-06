#!/usr/bin/env python3
"""Ionization-method card with an emulated signal on the published EXO-18-002 background model.
usage: make_card.py emu.npz name sigma_pb outdir [trig_syst]
Signal normalized to sigma_pb x L, so the limit on r is sigma95 / sigma_pb."""
import math
import os
import sys
import numpy as np
import ROOT
import summarize as S

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hepdata", "original", "ioniz", "AllCardsOnly") + "/"
REF = "Signal_tauPrime1e-1000"
LUMI = 101.0e3  # pb^-1
SIG_SYST = 0.15  # emulation closure on tau' 1e
NMAX = 1000.  # signal yield cap: Combine expected limits fail at r << 1e-2


def main(fn, name, sigma, outdir, trig_syst=None):
    os.makedirs(outdir, exist_ok=True)
    nom = S.bins(fn)
    tot = sum(v.sum() for v in nom.values())
    rscale = min(1., NMAX / (tot * sigma * LUMI)) if tot > 0 else 1.  # r_true = r_card * rscale
    sigma0, sigma = sigma, sigma * rscale
    if trig_syst is None:
        up, dn = S.bins(fn, "trig_up"), S.bins(fn, "trig_dn")
        tn = sum(v.sum() for v in nom.values())
        tu = sum(v.sum() for v in up.values()) / tn if tn > 0 else 1.
        td = sum(v.sum() for v in dn.values()) / tn if tn > 0 else 1.
    else:
        tu, td = 1 + trig_syst, 1 - trig_syst
    f = ROOT.TFile(D + "base.root")
    w = f.Get("w")
    for reg, eff in nom.items():
        ref = w.data("%s_%s" % (REF, reg))
        h = ref.Clone("Signal_%s_%s" % (name, reg))
        h.SetName("Signal_%s_%s" % (name, reg))
        for i in range(h.numEntries()):
            h.get(i)
            y = max(eff[i] * sigma * LUMI, 1e-9)
            h.set(i, y, 0.)
        getattr(w, "import")(h)
    fo = ROOT.TFile(os.path.join(outdir, "ws_%s.root" % name), "RECREATE")
    w.Write()
    fo.Close()
    card = []
    for l in open(D + "%s-0x0_area_card.txt" % REF):
        if l.split() and l.split()[1:2] == ["shapes"] and l.split()[0] != "lumi":
            continue  # signal shape systematics replaced by lnN below
        if l.startswith("kmax"):
            l = "kmax *\n"
        card.append(l.replace(REF, "Signal_" + name).replace("base.root", "ws_%s.root" % name))
    k = [i for i, l in enumerate(card) if l.startswith("lumi")][0]
    card.insert(k + 1, "sigeff lnN " + " ".join(["%.3f" % (1 + SIG_SYST)] * 4 + ["-"] * 4) + "\n")
    card.insert(k + 2, "trig lnN " + " ".join(["%.3f/%.3f" % (td, tu)] * 4 + ["-"] * 4) + "\n")
    ne = S.neff(fn)
    mc = 1. + min(1. / math.sqrt(ne), 1.) if ne > 0 else 2.  # MC statistics of the emulated signal
    card.insert(k + 3, "mcstat lnN " + " ".join(["%.3f" % mc] * 4 + ["-"] * 4) + "\n")
    open(os.path.join(outdir, "card_%s.txt" % name), "w").writelines(card)
    print("make_card: %s eps %.4f  N(sigma) %.3f  trig %.3f/%.3f  Neff %.1f  rscale %.6g" % (name, tot, tot * sigma0 * LUMI, td, tu, ne, rscale))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4], float(sys.argv[5]) if len(sys.argv) > 5 else None)
