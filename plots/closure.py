#!/usr/bin/env python3
"""tau' |Q| = 1e closure figures: sigma95 (emulated vs published), absolute efficiency and region fractions."""
import json
import math
import os
import re
import sys
import numpy as np
import ROOT
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "recast"))
import pstyle
import summarize
import qreco

OUT = os.environ.get("QP_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))  # outputs
BASE = OUT + "/"
PUB = os.path.join(qreco.HEPDATA, "ioniz", "AllCardsOnly")
LUMI = 101.0e3
cy = json.load(open(os.path.join(qreco.HERE, "card_yields.json")))
MS = [800, 1000, 1400, 1800, 2200]  # 2.6 TeV dropped: flat pedestal in the HEPData template


def lim(fn):
    text = open(fn).read()
    get = lambda key: float(re.search(key + r": r < (\S+)", text).group(1))
    return get("Observed Limit"), get("Expected 50.0%"), get("Expected 16.0%"), get("Expected 84.0%")


def mcstat(fn):
    """Emulation MC-statistics errors (weighted binomial) on the total efficiency, FAIL and last-bin fractions."""
    npz, col = summarize.load(fn)
    arr = npz["data"]
    ng = float(npz["ngen"])
    wgt = arr[:, col["trig"]] * ((arr[:, col["depth"]] >= summarize.NSTEP) & (arr[:, col["c_pt"]] > 200.))
    fail = arr[:, col["c_fpix"]] <= 0.9
    last = (~fail) & (arr[:, col["c_gstrip"]] > 0.4)
    wsum = wgt.sum()
    ferr = lambda mask: math.sqrt(np.sum(wgt ** 2 * (mask - (wgt * mask).sum() / wsum) ** 2)) / wsum
    return dict(eps_err=math.sqrt(np.sum(wgt ** 2)) / ng, fail_err=ferr(fail), last_err=ferr(last),
                rel_last=math.sqrt(np.sum((wgt * last) ** 2)) / (wgt * last).sum())


rows = []
for mass in MS:
    sr = cy[str(mass)]["sigma_ref"]
    po, pe, pm1, pp1 = lim("%s/pub_tp1e_%d.log" % (PUB, mass))
    eo, ee = lim("%s/cards/closure/lim_taup%d.log" % (BASE, mass))[:2]
    nexp = sr * LUMI
    pap = {reg: np.array(vals) / nexp for reg, vals in cy[str(mass)]["yields"].items()}
    emu = summarize.bins("%s/emu/taup_m%d_L1.npz" % (BASE, mass))
    tp, te = sum(vals.sum() for vals in pap.values()), sum(vals.sum() for vals in emu.values())
    st = mcstat("%s/emu/taup_m%d_L1.npz" % (BASE, mass))
    rows.append(dict(m=mass, T_pap=sum(sum(vals) for vals in cy[str(mass)]["yields"].values()), pub_obs=po * sr, pub_exp=pe * sr, pub_m1=pm1 * sr, pub_p1=pp1 * sr, **st, emu_obs=eo * sr, emu_exp=ee * sr, eps_pap=tp, eps_emu=te,
                     fail_pap=(pap["fail_LOW"].sum() + pap["fail_SIG"].sum()) / tp,
                     fail_emu=(emu["fail_LOW"].sum() + emu["fail_SIG"].sum()) / te,
                     last_pap=pap["pass_SIG"][-1] / tp, last_emu=emu["pass_SIG"][-1] / te))
json.dump(rows, open(os.path.join(pstyle.OUT, "closure_taup.json"), "w"), indent=1)
for row in rows:
    print("m %4d  sigma95 obs pub %.3g emu %.3g (%.2f)  exp pub %.3g emu %.3g (%.2f)  eps %.3f/%.3f  fail %.2f/%.2f  last %.2f/%.2f" % (
        row["m"], row["pub_obs"], row["emu_obs"], row["emu_obs"] / row["pub_obs"],
        row["pub_exp"], row["emu_exp"], row["emu_exp"] / row["pub_exp"],
        row["eps_pap"], row["eps_emu"], row["fail_pap"], row["fail_emu"], row["last_pap"], row["last_emu"]))

mtev = [row["m"] / 1000. for row in rows]
zeros = [0.] * len(rows)
pstyle.label("13 TeV, 101 fb^{-1}")
canv = pstyle.dicanvas("c_s95", 0.6, 2.4, 4e-5, 6e-4, 0.0, 2.2, "m_{#tau'} [TeV]", "#sigma_{95} [pb]", "Emul. / CMS", logy=True)
# emulation: MC-statistics band (relative error of the most sensitive bin)
for key, sty in (("emu_obs", 1), ("emu_exp", 2)):
    gband = pstyle.graph(mtev, [row[key] for row in rows], [row[key] * row["rel_last"] for row in rows])
    gband.SetFillColorAlpha(ROOT.kAzure + 1, 0.3)
    gband.Draw("3 same")
    ge_ = pstyle.graph(mtev, [row[key] for row in rows])
    ge_.SetLineColor(ROOT.kAzure + 1)
    ge_.SetLineStyle(sty)
    ge_.SetLineWidth(2)
    ge_.Draw("L same")
    if key == "emu_obs":
        ge, band = ge_, gband
    else:
        gee = ge_
gp = pstyle.graph(mtev, [row["pub_obs"] for row in rows])
gp.SetMarkerStyle(20)
gp.Draw("P same")
gpe = ROOT.TGraphAsymmErrors(len(rows))  # expected with its 68% band
for idx, row in enumerate(rows):
    gpe.SetPoint(idx, mtev[idx], row["pub_exp"])
    gpe.SetPointError(idx, 0, 0, row["pub_exp"] - row["pub_m1"], row["pub_p1"] - row["pub_exp"])
gpe.SetMarkerStyle(24)
gpe.Draw("P same")
pstyle.keep.append(gpe)
lg = pstyle.legend(0.40, 0.64, 0.93, 0.88, 0.030)
lg.AddEntry(gp, "CMS observed (HEPData)", "p")
lg.AddEntry(gpe, "CMS likelihood, expected (68% band)", "pe")
lg.AddEntry(ge, "Emulation, observed", "l")
lg.AddEntry(gee, "Emulation, expected", "l")
lg.AddEntry(band, "Emulation MC statistics", "f")
lg.Draw()
canv.cd(2)
gr = pstyle.graph(mtev, [row["emu_obs"] / row["pub_obs"] for row in rows],
                  [row["emu_obs"] / row["pub_obs"] * row["rel_last"] for row in rows])
gr.SetMarkerStyle(20)
gr.Draw("P same")
gre = ROOT.TGraphAsymmErrors(len(rows))  # MC statistics and the 68% band of the CMS expected, in quadrature
for idx, row in enumerate(rows):
    ratio = row["emu_exp"] / row["pub_exp"]
    up = ratio * math.hypot(row["rel_last"], (row["pub_exp"] - row["pub_m1"]) / row["pub_exp"])  # lower CMS limit raises the ratio
    dn = ratio * math.hypot(row["rel_last"], (row["pub_p1"] - row["pub_exp"]) / row["pub_exp"])
    gre.SetPoint(idx, mtev[idx] + 0.03, ratio)
    gre.SetPointError(idx, 0, 0, dn, up)
gre.SetMarkerStyle(24)
gre.Draw("P same")
pstyle.keep.append(gre)
ln = ROOT.TLine(0.6, 1., 2.4, 1.)
ln.SetLineStyle(2)
ln.Draw()
pstyle.keep.append(ln)
pstyle.save(canv, "closure_sigma95")

c2 = pstyle.dicanvas("c_frac", 0.6, 2.4, 0., 1.45, 0.3, 2.4, "m_{#tau'} [TeV]", "Efficiency or fraction", "CMS / Emul.")
items = [("last", "Fraction in PASS, G_{i}^{Strips} > 0.4", ROOT.kAzure + 1), ("eps", "Total efficiency", ROOT.kBlack),
         ("fail", "Fraction in FAIL", ROOT.kOrange + 7)]
lg2 = pstyle.legend(0.24, 0.72, 0.93, 0.88, 0.032)  # above the data
for qty, title, col in items:
    gband = pstyle.graph(mtev, [row[qty + "_emu"] for row in rows], [row[qty + "_err"] for row in rows])
    gband.SetFillColorAlpha(col, 0.3)
    gband.Draw("3 same")
    ge = pstyle.graph(mtev, [row[qty + "_emu"] for row in rows])
    ge.SetLineColor(col)
    ge.SetLineWidth(2)
    ge.Draw("L same")
    gp = pstyle.graph(mtev, [row[qty + "_pap"] for row in rows])
    gp.SetMarkerStyle(20)
    gp.SetMarkerColor(col)
    gp.Draw("P same")
    lg2.AddEntry(ge, title, "l")
lg2.Draw()
pstyle.note("Points: CMS likelihood; lines: emulation (band: MC stat.)", 0.26, 0.68, 0.030)
c2.cd(2)
sb = ROOT.TBox(0.6, 0.85, 2.4, 1.15)  # assigned signal-efficiency uncertainty
sb.SetFillColorAlpha(ROOT.kGray, 0.6)
sb.Draw()
pstyle.keep.append(sb)
for idx, (qty, title, col) in enumerate(items):
    ratios = [row[qty + "_pap"] / row[qty + "_emu"] for row in rows]
    gratio = pstyle.graph([xv + 0.015 * (idx - 1) for xv in mtev], ratios,
                          [qi * row[qty + "_err"] / row[qty + "_emu"] for qi, row in zip(ratios, rows)])
    gratio.SetMarkerStyle(20)
    gratio.SetMarkerColor(col)
    gratio.SetLineColor(col)
    gratio.Draw("P same")
ln2 = ROOT.TLine(0.6, 1., 2.4, 1.)
ln2.SetLineStyle(2)
ln2.Draw()
pstyle.keep.append(ln2)
pstyle.note("gray: #pm15% assigned uncertainty", 0.20, 0.80, 0.07)
pstyle.save(c2, "closure_fractions")

# variant: Poisson errors on the published yields (expected signal events at sigma_ref x L), ratio CMS / emulation
c3 = pstyle.dicanvas("c_frac_poisson", 0.6, 2.4, 0., 1.45, 0., 3., "m_{#tau'} [TeV]", "Efficiency or fraction", "CMS / Emul.")
lg3 = pstyle.legend(0.24, 0.72, 0.93, 0.88, 0.032)
pts = {}
for idx, (qty, title, col) in enumerate(items):
    gband = pstyle.graph(mtev, [row[qty + "_emu"] for row in rows], [row[qty + "_err"] for row in rows])
    gband.SetFillColorAlpha(col, 0.3)
    gband.Draw("3 same")
    ge = pstyle.graph(mtev, [row[qty + "_emu"] for row in rows])
    ge.SetLineColor(col)
    ge.SetLineWidth(2)
    ge.Draw("L same")
    perr = []
    for row in rows:
        frac, total = row[qty + "_pap"], row["T_pap"]
        perr.append(frac / math.sqrt(total) if qty == "eps" else math.sqrt(frac * (1 - frac) / total))  # Poisson on the card yields
    gp = pstyle.graph([xv + 0.012 * (idx - 1) for xv in mtev], [row[qty + "_pap"] for row in rows], perr)
    gp.SetMarkerStyle(20)
    gp.SetMarkerColor(col)
    gp.SetLineColor(col)
    gp.Draw("P same")
    pts[qty] = (perr, col)
    lg3.AddEntry(ge, title, "l")
lg3.Draw()
pstyle.note("Points: CMS likelihood, #sqrt{N} on the expected yields", 0.26, 0.68, 0.028)
c3.cd(2)
for idx, (qty, title, col) in enumerate(items):
    perr, col = pts[qty]
    ratios = [row[qty + "_pap"] / row[qty + "_emu"] for row in rows]
    eq = [qi * math.hypot(pe / row[qty + "_pap"], row[qty + "_err"] / row[qty + "_emu"]) for qi, pe, row in zip(ratios, perr, rows)]
    gratio = pstyle.graph([xv + 0.012 * (idx - 1) for xv in mtev], ratios, eq)
    gratio.SetMarkerStyle(20)
    gratio.SetMarkerColor(col)
    gratio.SetLineColor(col)
    gratio.Draw("P same")
ln3 = ROOT.TLine(0.6, 1., 2.4, 1.)
ln3.SetLineStyle(2)
ln3.Draw()
pstyle.keep.append(ln3)
pstyle.save(c3, "closure_fractions_poisson")
