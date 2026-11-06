#!/usr/bin/env python3
"""Per-bin signal yields of the published tau' 1e ionization cards + sigma_ref = sigma95(Fig. 7) / r_obs."""
import json, os, re, yaml, ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hepdata", "original") + "/"
f7 = yaml.safe_load(open(D + "figure_7_left.yaml"))
ms = [round(1000 * float(v["value"])) for v in f7["independent_variables"][0]["values"]]
s95 = dict(zip(ms, [float(v["value"]) for v in f7["dependent_variables"][1]["values"]]))
out = {}
for m in (200, 800, 1000, 1400, 1800, 2200, 2600):
    card = D + "ioniz/AllCardsOnly/Signal_tauPrime1e-%d-0x0_area_card.txt" % m
    fn = [l.split()[3] for l in open(card) if l.startswith("shapes") and "Signal" in l][0]
    w = ROOT.TFile(D + "ioniz/AllCardsOnly/" + fn).Get("w")
    r = float(re.search(r"Observed Limit: r < (\S+)", open(D + "ioniz/AllCardsOnly/pub_tp1e_%d.log" % m).read()).group(1))
    y = {}
    for reg in ("fail_LOW", "fail_SIG", "pass_LOW", "pass_SIG"):
        h = w.data("Signal_tauPrime1e-%d_%s" % (m, reg))
        y[reg] = [h.weight(i) if (h.get(i), True)[1] else 0 for i in range(h.numEntries())]
        y[reg] = []
        for i in range(h.numEntries()):
            h.get(i); y[reg].append(h.weight())
    out[m] = dict(r_obs=r, s95=s95.get(m), sigma_ref=s95[m] / r if m in s95 else None, yields=y,
                  total=sum(sum(v) for v in y.values()))
    print(m, "r_obs %.4g" % r, "sigma_ref", out[m]["sigma_ref"], "total %.3f" % out[m]["total"], "pass_SIG", [round(v, 3) for v in y["pass_SIG"]])
json.dump(out, open("card_yields.json", "w"), indent=1)
