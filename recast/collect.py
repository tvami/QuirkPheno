#!/usr/bin/env python3
"""Collect quirk limits: r = sigma95 / sigma_LO (observed, expected, +-1 sigma) per method and the best expected."""
import glob, json, os, re, sys
OUT = os.environ.get("QP_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))  # outputs
CARDDIR = os.path.join(OUT, "cards", sys.argv[1] if len(sys.argv) > 1 else "nom")
FAM = sys.argv[2] if len(sys.argv) > 2 else "doublet"
NEFF_MIN = 3.  # minimum effective MC events in the signal region
MODELS = {"doublet": ("quirkb", "quirk"), "singlet": ("quirksb", "quirks")}[FAM]  # bound-filtered sample first


def lim(fn):
    if not os.path.exists(fn):
        return None
    text = open(fn).read()
    try:
        grab = lambda label: float(re.search(label + r": r < (\S+)", text).group(1))
        return dict(obs=grab("Observed Limit"), exp=grab("Expected 50.0%"), m1=grab("Expected 16.0%"), p1=grab("Expected 84.0%"))
    except AttributeError:
        return None


res, seen = [], set()
# the bound-pair samples (quirkb, GEN-filtered, high statistics) replace quirk where both exist
for xs in sorted(glob.glob(CARDDIR + "/xs_%s_*.txt" % MODELS[0])) + sorted(glob.glob(CARDDIR + "/xs_%s_*.txt" % MODELS[1])):
    mass, lam, sigma = open(xs).read().split()
    name = os.path.basename(xs)[3:-4]
    if (mass, lam) in seen:
        continue
    seen.add((mass, lam))
    ion, mas = lim("%s/lim_%s.log" % (CARDDIR, name)), lim("%s/lim_mass_%s.log" % (CARDDIR, name))
    best = ion  # ionization method only (mass method dropped, 2026-11-02)
    mk = "%s/card_%s.mklog" % (CARDDIR, name)
    nm = re.findall(r"Neff ([\d.]+)", open(mk).read()) if os.path.exists(mk) else []
    if nm and float(nm[-1]) < NEFF_MIN:  # too few MC events to trust a limit
        best = None
    rs = re.findall(r"rscale ([\d.e+-]+)", open(mk).read()) if os.path.exists(mk) else []
    if rs and ion:  # card built at a capped signal yield
        ion = {key: val * float(rs[-1]) for key, val in ion.items()}
        best = best and ion
    res.append(dict(m=int(mass), L=float(lam), xsec=float(sigma), ion=ion, mass=mas, best=best, model=name.split("_")[0],
                    method="ion"))
json.dump(res, open(CARDDIR + ("/limits.json" if FAM == "doublet" else "/limits_%s.json" % FAM), "w"), indent=1)
for row in sorted(res, key=lambda row: (row["L"], row["m"])):
    fmt = lambda entry: "%8.3g" % entry["obs"] if entry else "       -"
    print("L %7g m %5d sigma %.3g  r_obs ion %s mass %s  best(%s) obs %s exp %s" % (row["L"], row["m"], row["xsec"],
          fmt(row["ion"]), fmt(row["mass"]), row["method"], fmt(row["best"]),
          "%8.3g" % row["best"]["exp"] if row["best"] else "-"))
