#!/usr/bin/env python3
"""Per-bin signal yields of the published tau' 1e ionization cards + sigma_ref = sigma95(Fig. 7) / r_obs."""
import json, os, re, yaml, ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
HEPDATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hepdata", "original") + "/"
CARDS = HEPDATA + "ioniz/AllCardsOnly/"
fig7 = yaml.safe_load(open(HEPDATA + "figure_7_left.yaml"))
fig7_masses = [round(1000 * float(value["value"])) for value in fig7["independent_variables"][0]["values"]]
sigma95 = dict(zip(fig7_masses, [float(value["value"]) for value in fig7["dependent_variables"][1]["values"]]))
out = {}
for mass in (200, 800, 1000, 1400, 1800, 2200):
    card = CARDS + "Signal_tauPrime1e-%d-0x0_area_card.txt" % mass
    ws_file = [line.split()[3] for line in open(card) if line.startswith("shapes") and "Signal" in line][0]
    workspace = ROOT.TFile(CARDS + ws_file).Get("w")
    r_obs = float(re.search(r"Observed Limit: r < (\S+)", open(CARDS + "pub_tp1e_%d.log" % mass).read()).group(1))
    yields = {}
    for region in ("fail_LOW", "fail_SIG", "pass_LOW", "pass_SIG"):
        data = workspace.data("Signal_tauPrime1e-%d_%s" % (mass, region))
        yields[region] = []
        for index in range(data.numEntries()):
            data.get(index)
            yields[region].append(data.weight())
    out[mass] = dict(r_obs=r_obs, s95=sigma95.get(mass), sigma_ref=sigma95[mass] / r_obs if mass in sigma95 else None, yields=yields,
                     total=sum(sum(values) for values in yields.values()))
    print(mass, "r_obs %.4g" % r_obs, "sigma_ref", out[mass]["sigma_ref"], "total %.3f" % out[mass]["total"],
          "pass_SIG", [round(value, 3) for value in yields["pass_SIG"]])
json.dump(out, open("card_yields.json", "w"), indent=1)
