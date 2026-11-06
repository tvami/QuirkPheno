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
import pstyle as P
import summarize as S

B = "/ceph/cms/store/user/tvami/HSCP-Quirks/QuirkPheno/"
PUB = os.path.join(S.Q.HEPDATA, "ioniz", "AllCardsOnly")
L = 101.0e3
cy = json.load(open(os.path.join(S.Q.HERE, "card_yields.json")))
MS = [800, 1000, 1400, 1800, 2200]  # 2.6 TeV dropped: flat pedestal in the HEPData template


def lim(fn):
    t = open(fn).read()
    g = lambda k: float(re.search(k + r": r < (\S+)", t).group(1))
    return g("Observed Limit"), g("Expected 50.0%"), g("Expected 16.0%"), g("Expected 84.0%")


def mcstat(fn):
    """Emulation MC-statistics errors (weighted binomial) on the total efficiency, FAIL and last-bin fractions."""
    d, k = S.load(fn)
    x = d["data"]
    ng = float(d["ngen"])
    w = x[:, k["trig"]] * ((x[:, k["depth"]] >= S.NSTEP) & (x[:, k["c_pt"]] > 200.))
    fail = x[:, k["c_fpix"]] <= 0.9
    last = (~fail) & (x[:, k["c_gstrip"]] > 0.4)
    T = w.sum()
    ferr = lambda a: math.sqrt(np.sum(w ** 2 * (a - (w * a).sum() / T) ** 2)) / T
    return dict(eps_err=math.sqrt(np.sum(w ** 2)) / ng, fail_err=ferr(fail), last_err=ferr(last),
                rel_last=math.sqrt(np.sum((w * last) ** 2)) / (w * last).sum())


rows = []
for m in MS:
    sr = cy[str(m)]["sigma_ref"]
    po, pe, pm1, pp1 = lim("%s/pub_tp1e_%d.log" % (PUB, m))
    eo, ee = lim("%s/cards/closure/lim_taup%d.log" % (B, m))[:2]
    n = sr * L
    pap = {k: np.array(v) / n for k, v in cy[str(m)]["yields"].items()}
    emu = S.bins("%s/emu/taup_m%d_L1.npz" % (B, m))
    tp, te = sum(v.sum() for v in pap.values()), sum(v.sum() for v in emu.values())
    st = mcstat("%s/emu/taup_m%d_L1.npz" % (B, m))
    rows.append(dict(m=m, T_pap=sum(sum(v) for v in cy[str(m)]["yields"].values()), pub_obs=po * sr, pub_exp=pe * sr, pub_m1=pm1 * sr, pub_p1=pp1 * sr, **st, emu_obs=eo * sr, emu_exp=ee * sr, eps_pap=tp, eps_emu=te,
                     fail_pap=(pap["fail_LOW"].sum() + pap["fail_SIG"].sum()) / tp,
                     fail_emu=(emu["fail_LOW"].sum() + emu["fail_SIG"].sum()) / te,
                     last_pap=pap["pass_SIG"][-1] / tp, last_emu=emu["pass_SIG"][-1] / te))
json.dump(rows, open(os.path.join(P.OUT, "closure_taup.json"), "w"), indent=1)
for r in rows:
    print("m %4d  sigma95 obs pub %.3g emu %.3g (%.2f)  exp pub %.3g emu %.3g (%.2f)  eps %.3f/%.3f  fail %.2f/%.2f  last %.2f/%.2f" % (
        r["m"], r["pub_obs"], r["emu_obs"], r["emu_obs"] / r["pub_obs"], r["pub_exp"], r["emu_exp"], r["emu_exp"] / r["pub_exp"],
        r["eps_pap"], r["eps_emu"], r["fail_pap"], r["fail_emu"], r["last_pap"], r["last_emu"]))

x = [r["m"] / 1000. for r in rows]
z = [0.] * len(rows)
P.label("13 TeV, 101 fb^{-1}")
c = P.dicanvas("c_s95", 0.6, 2.4, 4e-5, 6e-4, 0.0, 2.2, "m_{#tau'} [TeV]", "#sigma_{95} [pb]", "Emul. / CMS", logy=True)
# emulation: MC-statistics band (relative error of the most sensitive bin)
for key, sty in (("emu_obs", 1), ("emu_exp", 2)):
    b = P.graph(x, [r[key] for r in rows], [r[key] * r["rel_last"] for r in rows])
    b.SetFillColorAlpha(ROOT.kAzure + 1, 0.3)
    b.Draw("3 same")
    ge_ = P.graph(x, [r[key] for r in rows])
    ge_.SetLineColor(ROOT.kAzure + 1)
    ge_.SetLineStyle(sty)
    ge_.SetLineWidth(2)
    ge_.Draw("L same")
    if key == "emu_obs":
        ge, band = ge_, b
    else:
        gee = ge_
gp = P.graph(x, [r["pub_obs"] for r in rows])
gp.SetMarkerStyle(20)
gp.Draw("P same")
gpe = ROOT.TGraphAsymmErrors(len(rows))  # expected with its 68% band
for i, r in enumerate(rows):
    gpe.SetPoint(i, x[i], r["pub_exp"])
    gpe.SetPointError(i, 0, 0, r["pub_exp"] - r["pub_m1"], r["pub_p1"] - r["pub_exp"])
gpe.SetMarkerStyle(24)
gpe.Draw("P same")
P.keep.append(gpe)
lg = P.legend(0.40, 0.64, 0.93, 0.88, 0.030)
lg.AddEntry(gp, "CMS observed (HEPData)", "p")
lg.AddEntry(gpe, "CMS likelihood, expected (68% band)", "pe")
lg.AddEntry(ge, "Emulation, observed", "l")
lg.AddEntry(gee, "Emulation, expected", "l")
lg.AddEntry(band, "Emulation MC statistics", "f")
lg.Draw()
c.cd(2)
gr = P.graph(x, [r["emu_obs"] / r["pub_obs"] for r in rows], [r["emu_obs"] / r["pub_obs"] * r["rel_last"] for r in rows])
gr.SetMarkerStyle(20)
gr.Draw("P same")
gre = ROOT.TGraphAsymmErrors(len(rows))  # MC statistics and the 68% band of the CMS expected, in quadrature
for i, r in enumerate(rows):
    q = r["emu_exp"] / r["pub_exp"]
    up = q * math.hypot(r["rel_last"], (r["pub_exp"] - r["pub_m1"]) / r["pub_exp"])  # lower CMS limit raises the ratio
    dn = q * math.hypot(r["rel_last"], (r["pub_p1"] - r["pub_exp"]) / r["pub_exp"])
    gre.SetPoint(i, x[i] + 0.03, q)
    gre.SetPointError(i, 0, 0, dn, up)
gre.SetMarkerStyle(24)
gre.Draw("P same")
P.keep.append(gre)
ln = ROOT.TLine(0.6, 1., 2.4, 1.)
ln.SetLineStyle(2)
ln.Draw()
P.keep.append(ln)
P.save(c, "closure_sigma95")

c2 = P.dicanvas("c_frac", 0.6, 2.4, 0., 1.45, 0.3, 2.4, "m_{#tau'} [TeV]", "Efficiency or fraction", "CMS / Emul.")
items = [("last", "Fraction in PASS, G_{i}^{Strips} > 0.4", ROOT.kAzure + 1), ("eps", "Total efficiency", ROOT.kBlack),
         ("fail", "Fraction in FAIL", ROOT.kOrange + 7)]
lg2 = P.legend(0.24, 0.72, 0.93, 0.88, 0.032)  # above the data
for k, t, col in items:
    b = P.graph(x, [r[k + "_emu"] for r in rows], [r[k + "_err"] for r in rows])
    b.SetFillColorAlpha(col, 0.3)
    b.Draw("3 same")
    ge = P.graph(x, [r[k + "_emu"] for r in rows])
    ge.SetLineColor(col)
    ge.SetLineWidth(2)
    ge.Draw("L same")
    gp = P.graph(x, [r[k + "_pap"] for r in rows])
    gp.SetMarkerStyle(20)
    gp.SetMarkerColor(col)
    gp.Draw("P same")
    lg2.AddEntry(ge, t, "l")
lg2.Draw()
P.note("Points: CMS likelihood; lines: emulation (band: MC stat.)", 0.26, 0.68, 0.030)
c2.cd(2)
sb = ROOT.TBox(0.6, 0.85, 2.4, 1.15)  # assigned signal-efficiency uncertainty
sb.SetFillColorAlpha(ROOT.kGray, 0.6)
sb.Draw()
P.keep.append(sb)
for j, (k, t, col) in enumerate(items):
    q = [r[k + "_pap"] / r[k + "_emu"] for r in rows]
    g = P.graph([v + 0.015 * (j - 1) for v in x], q, [qi * r[k + "_err"] / r[k + "_emu"] for qi, r in zip(q, rows)])
    g.SetMarkerStyle(20)
    g.SetMarkerColor(col)
    g.SetLineColor(col)
    g.Draw("P same")
ln2 = ROOT.TLine(0.6, 1., 2.4, 1.)
ln2.SetLineStyle(2)
ln2.Draw()
P.keep.append(ln2)
P.note("gray: #pm15% assigned uncertainty", 0.20, 0.80, 0.07)
P.save(c2, "closure_fractions")

# variant: Poisson errors on the published yields (expected signal events at sigma_ref x L), ratio CMS / emulation
c3 = P.dicanvas("c_frac_poisson", 0.6, 2.4, 0., 1.45, 0., 3., "m_{#tau'} [TeV]", "Efficiency or fraction", "CMS / Emul.")
lg3 = P.legend(0.24, 0.72, 0.93, 0.88, 0.032)
pts = {}
for j, (k, t, col) in enumerate(items):
    b = P.graph(x, [r[k + "_emu"] for r in rows], [r[k + "_err"] for r in rows])
    b.SetFillColorAlpha(col, 0.3)
    b.Draw("3 same")
    ge = P.graph(x, [r[k + "_emu"] for r in rows])
    ge.SetLineColor(col)
    ge.SetLineWidth(2)
    ge.Draw("L same")
    perr = []
    for r in rows:
        f, T = r[k + "_pap"], r["T_pap"]
        perr.append(f / math.sqrt(T) if k == "eps" else math.sqrt(f * (1 - f) / T))  # Poisson on the card yields
    gp = P.graph([v + 0.012 * (j - 1) for v in x], [r[k + "_pap"] for r in rows], perr)
    gp.SetMarkerStyle(20)
    gp.SetMarkerColor(col)
    gp.SetLineColor(col)
    gp.Draw("P same")
    pts[k] = (perr, col)
    lg3.AddEntry(ge, t, "l")
lg3.Draw()
P.note("Points: CMS likelihood, #sqrt{N} on the expected yields", 0.26, 0.68, 0.028)
c3.cd(2)
for j, (k, t, col) in enumerate(items):
    perr, col = pts[k]
    q = [r[k + "_pap"] / r[k + "_emu"] for r in rows]
    eq = [qi * math.hypot(pe / r[k + "_pap"], r[k + "_err"] / r[k + "_emu"]) for qi, pe, r in zip(q, perr, rows)]
    g = P.graph([v + 0.012 * (j - 1) for v in x], q, eq)
    g.SetMarkerStyle(20)
    g.SetMarkerColor(col)
    g.SetLineColor(col)
    g.Draw("P same")
ln3 = ROOT.TLine(0.6, 1., 2.4, 1.)
ln3.SetLineStyle(2)
ln3.Draw()
P.keep.append(ln3)
P.save(c3, "closure_fractions_poisson")
