#!/usr/bin/env python3
"""Cutflow and ionization-method bin efficiencies from emulate.py output. usage: summarize.py emu.npz [paper_mass]"""
import sys
import numpy as np
import yaml
import qreco

LOW = [0., 0.05, 0.1, 0.15, 0.2, 0.25]
SIG = [0.25, 0.3, 0.35, 0.4, 1.0001]
NSTEP = 14  # candidate steps after the trigger


def load(fn):
    npz = np.load(fn)
    col = {name: idx for idx, name in enumerate(npz["keys"])}
    return npz, col


def cutflow(fn, trig="trig"):
    npz, col = load(fn)
    data = npz["data"]
    weight = data[:, col[trig]]
    ngen = float(npz["ngen"])
    out = [1., weight.sum() / ngen]
    for step in range(1, NSTEP + 1):
        out.append(weight[data[:, col["depth"]] >= step].sum() / ngen)
    return np.array(out)


def bins(fn, trig="trig", ptmin=200.):
    """Efficiency per (fail/pass, LOW/SIG) Gstrip bin, as in the datacards (pT > 200)."""
    npz, col = load(fn)
    data = npz["data"]
    sel = (data[:, col["depth"]] >= NSTEP) & (data[:, col["c_pt"]] > ptmin)
    weight = data[:, col[trig]] * sel / float(npz["ngen"])
    gstrip, fpix = data[:, col["c_gstrip"]], data[:, col["c_fpix"]]
    res = {}
    for reg, cond in (("fail", fpix <= 0.9), ("pass", fpix > 0.9)):
        res[reg + "_LOW"] = np.histogram(gstrip, LOW, weights=weight * cond)[0]
        res[reg + "_SIG"] = np.histogram(gstrip, SIG, weights=weight * cond)[0]
    return res


def paper_cutflow(mass):
    doc = yaml.safe_load(open(qreco.HEPDATA + "/efficiency_tau1e.yaml"))
    for dv in doc["dependent_variables"]:
        if dv["header"]["name"].strip() == str(mass):
            return np.array([float(val["value"]) for val in dv["values"]])
    return None


if __name__ == "__main__":
    import emulate
    fn = sys.argv[1]
    cf = cutflow(fn)
    pm = paper_cutflow(int(sys.argv[2])) if len(sys.argv) > 2 else None
    for idx, cut in enumerate(emulate.CUTS):
        print("%-14s %.3f %s" % (cut, cf[idx], "  paper %.3f" % pm[idx] if pm is not None else ""))
    for reg, eff in bins(fn).items():
        print(reg, np.round(eff, 4), "sum %.4f" % eff.sum())


def neff(fn, trig="trig", ptmin=200.):
    """Effective number of MC events in the signal region, (sum w)^2 / sum w^2."""
    npz, col = load(fn)
    data = npz["data"]
    weight = data[:, col[trig]] * ((data[:, col["depth"]] >= NSTEP) & (data[:, col["c_pt"]] > ptmin))
    s2 = float((weight ** 2).sum())
    return float(weight.sum()) ** 2 / s2 if s2 > 0 else 0.
