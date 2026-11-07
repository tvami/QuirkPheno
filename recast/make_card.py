#!/usr/bin/env python3
"""Ionization-method card with an emulated signal on the published EXO-18-002 background model.
usage: make_card.py emu.npz name sigma_pb outdir [trig_syst]
Signal normalized to sigma_pb x L, so the limit on r is sigma95 / sigma_pb."""
import math
import os
import sys
import numpy as np
import ROOT
import summarize

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
CARDS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hepdata", "original", "ioniz", "AllCardsOnly") + "/"
REF = "Signal_tauPrime1e-1000"
LUMI = 101.0e3  # pb^-1
SIG_SYST = 0.15  # emulation closure on tau' 1e
NMAX = 1000.  # signal yield cap: Combine expected limits fail at r << 1e-2


def main(fn, name, sigma, outdir, trig_syst=None):
    os.makedirs(outdir, exist_ok=True)
    nom = summarize.bins(fn)
    tot = sum(vals.sum() for vals in nom.values())
    rscale = min(1., NMAX / (tot * sigma * LUMI)) if tot > 0 else 1.  # r_true = r_card * rscale
    sigma0, sigma = sigma, sigma * rscale
    if trig_syst is None:
        up, dn = summarize.bins(fn, "trig_up"), summarize.bins(fn, "trig_dn")
        tn = sum(vals.sum() for vals in nom.values())
        tu = sum(vals.sum() for vals in up.values()) / tn if tn > 0 else 1.
        td = sum(vals.sum() for vals in dn.values()) / tn if tn > 0 else 1.
    else:
        tu, td = 1 + trig_syst, 1 - trig_syst
    tfile = ROOT.TFile(CARDS + "base.root")
    wspace = tfile.Get("w")
    for reg, eff in nom.items():
        ref = wspace.data("%s_%s" % (REF, reg))
        hist = ref.Clone("Signal_%s_%s" % (name, reg))
        hist.SetName("Signal_%s_%s" % (name, reg))
        for idx in range(hist.numEntries()):
            hist.get(idx)
            nsig = max(eff[idx] * sigma * LUMI, 1e-9)
            hist.set(idx, nsig, 0.)
        getattr(wspace, "import")(hist)
    fo = ROOT.TFile(os.path.join(outdir, "ws_%s.root" % name), "RECREATE")
    wspace.Write()
    fo.Close()
    card = []
    for line in open(CARDS + "%s-0x0_area_card.txt" % REF):
        if line.split() and line.split()[1:2] == ["shapes"] and line.split()[0] != "lumi":
            continue  # signal shape systematics replaced by lnN below
        if line.startswith("kmax"):
            line = "kmax *\n"
        card.append(line.replace(REF, "Signal_" + name).replace("base.root", "ws_%s.root" % name))
    pos = [idx for idx, line in enumerate(card) if line.startswith("lumi")][0]
    card.insert(pos + 1, "sigeff lnN " + " ".join(["%.3f" % (1 + SIG_SYST)] * 4 + ["-"] * 4) + "\n")
    card.insert(pos + 2, "trig lnN " + " ".join(["%.3f/%.3f" % (td, tu)] * 4 + ["-"] * 4) + "\n")
    ne = summarize.neff(fn)
    mc = 1. + min(1. / math.sqrt(ne), 1.) if ne > 0 else 2.  # MC statistics of the emulated signal
    card.insert(pos + 3, "mcstat lnN " + " ".join(["%.3f" % mc] * 4 + ["-"] * 4) + "\n")
    open(os.path.join(outdir, "card_%s.txt" % name), "w").writelines(card)
    print("make_card: %s eps %.4f  N(sigma) %.3f  trig %.3f/%.3f  Neff %.1f  rscale %.6g" % (name, tot, tot * sigma0 * LUMI, td, tu, ne, rscale))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4], float(sys.argv[5]) if len(sys.argv) > 5 else None)
