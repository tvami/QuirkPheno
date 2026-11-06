#!/usr/bin/env python3
"""Cutflow and ionization-method bin efficiencies from emulate.py output. usage: summarize.py emu.npz [paper_mass]"""
import sys
import numpy as np
import yaml
import qreco as Q

LOW = [0., 0.05, 0.1, 0.15, 0.2, 0.25]
SIG = [0.25, 0.3, 0.35, 0.4, 1.0001]
NSTEP = 14  # candidate steps after the trigger


def load(fn):
    d = np.load(fn)
    k = {s: i for i, s in enumerate(d["keys"])}
    return d, k


def cutflow(fn, trig="trig"):
    d, k = load(fn)
    x = d["data"]
    w = x[:, k[trig]]
    ngen = float(d["ngen"])
    out = [1., w.sum() / ngen]
    for j in range(1, NSTEP + 1):
        out.append(w[x[:, k["depth"]] >= j].sum() / ngen)
    return np.array(out)


def bins(fn, trig="trig", ptmin=200.):
    """Efficiency per (fail/pass, LOW/SIG) Gstrip bin, as in the datacards (pT > 200)."""
    d, k = load(fn)
    x = d["data"]
    sel = (x[:, k["depth"]] >= NSTEP) & (x[:, k["c_pt"]] > ptmin)
    w = x[:, k[trig]] * sel / float(d["ngen"])
    g, f = x[:, k["c_gstrip"]], x[:, k["c_fpix"]]
    res = {}
    for reg, cond in (("fail", f <= 0.9), ("pass", f > 0.9)):
        res[reg + "_LOW"] = np.histogram(g, LOW, weights=w * cond)[0]
        res[reg + "_SIG"] = np.histogram(g, SIG, weights=w * cond)[0]
    return res


def paper_cutflow(mass):
    d = yaml.safe_load(open(Q.HEPDATA + "/efficiency_tau1e.yaml"))
    for dv in d["dependent_variables"]:
        if dv["header"]["name"].strip() == str(mass):
            return np.array([float(v["value"]) for v in dv["values"]])
    return None


if __name__ == "__main__":
    import emulate
    fn = sys.argv[1]
    cf = cutflow(fn)
    pm = paper_cutflow(int(sys.argv[2])) if len(sys.argv) > 2 else None
    for i, c in enumerate(emulate.CUTS):
        print("%-14s %.3f %s" % (c, cf[i], "  paper %.3f" % pm[i] if pm is not None else ""))
    for r, v in bins(fn).items():
        print(r, np.round(v, 4), "sum %.4f" % v.sum())


def neff(fn, trig="trig", ptmin=200.):
    """Effective number of MC events in the signal region, (sum w)^2 / sum w^2."""
    d, k = load(fn)
    x = d["data"]
    w = x[:, k[trig]] * ((x[:, k["depth"]] >= NSTEP) & (x[:, k["c_pt"]] > ptmin))
    s2 = float((w ** 2).sum())
    return float(w.sum()) ** 2 / s2 if s2 > 0 else 0.
