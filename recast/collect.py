#!/usr/bin/env python3
"""Collect quirk limits: r = sigma95 / sigma_LO (observed, expected, +-1 sigma) per method and the best expected."""
import glob, json, os, re, sys
OUT = os.environ.get("QP_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))  # outputs
O = os.path.join(OUT, "cards", sys.argv[1] if len(sys.argv) > 1 else "nom")
FAM = sys.argv[2] if len(sys.argv) > 2 else "doublet"
NEFF_MIN = 3.  # minimum effective MC events in the signal region
MODELS = {"doublet": ("quirkb", "quirk"), "singlet": ("quirksb", "quirks")}[FAM]  # bound-filtered sample first


def lim(fn):
    if not os.path.exists(fn):
        return None
    t = open(fn).read()
    try:
        g = lambda k: float(re.search(k + r": r < (\S+)", t).group(1))
        return dict(obs=g("Observed Limit"), exp=g("Expected 50.0%"), m1=g("Expected 16.0%"), p1=g("Expected 84.0%"))
    except AttributeError:
        return None


res, seen = [], set()
# the bound-pair samples (quirkb, GEN-filtered, high statistics) replace quirk where both exist
for xs in sorted(glob.glob(O + "/xs_%s_*.txt" % MODELS[0])) + sorted(glob.glob(O + "/xs_%s_*.txt" % MODELS[1])):
    m, L, s = open(xs).read().split()
    n = os.path.basename(xs)[3:-4]
    if (m, L) in seen:
        continue
    seen.add((m, L))
    ion, mas = lim("%s/lim_%s.log" % (O, n)), lim("%s/lim_mass_%s.log" % (O, n))
    best = ion  # ionization method only (mass method dropped, 2026-11-02)
    mk = "%s/card_%s.mklog" % (O, n)
    nm = re.findall(r"Neff ([\d.]+)", open(mk).read()) if os.path.exists(mk) else []
    if nm and float(nm[-1]) < NEFF_MIN:  # too few MC events to trust a limit
        best = None
    rs = re.findall(r"rscale ([\d.e+-]+)", open(mk).read()) if os.path.exists(mk) else []
    if rs and ion:  # card built at a capped signal yield
        ion = {k: v * float(rs[-1]) for k, v in ion.items()}
        best = best and ion
    res.append(dict(m=int(m), L=float(L), xsec=float(s), ion=ion, mass=mas, best=best, model=n.split("_")[0],
                    method="ion"))
json.dump(res, open(O + ("/limits.json" if FAM == "doublet" else "/limits_%s.json" % FAM), "w"), indent=1)
for r in sorted(res, key=lambda r: (r["L"], r["m"])):
    f = lambda x: "%8.3g" % x["obs"] if x else "       -"
    print("L %7g m %5d sigma %.3g  r_obs ion %s mass %s  best(%s) obs %s exp %s" % (r["L"], r["m"], r["xsec"], f(r["ion"]), f(r["mass"]), r["method"], f(r["best"]), "%8.3g" % r["best"]["exp"] if r["best"] else "-"))
