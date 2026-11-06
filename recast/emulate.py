#!/usr/bin/env python3
"""EXO-18-002 selection emulated on quirksim output. usage: emulate.py <outdir> <tag-glob> <out.npz>
Writes the per-event cutflow and the best candidate (pt, p, Gstrip, Fpix, Ih) for every trigger variation."""
import glob
import math
import sys
import numpy as np
import os
import qreco as Q

MERGE = float(os.environ.get("QP_MERGE", 1.))  # two-hit resolution scale (systematic)

CUTS = ["All", "Trigger", "pT>55", "|eta|<1", "Npix>=2", "fvalid>0.8", "NdEdx>=10", "HighPurity", "chi2<5",
        "dz<0.1", "dxy<0.02", "MiniIso<0.02", "TkIso<15", "E/p<0.3", "sigpt<8e-4", "Fpix>0.3"]


def clusters(a, i, thick, mip, rng, beta_q):
    """Tracker clusters of one event: per (vol) merge the two quirks' hits within the two-hit resolution."""
    det = a["h_det"][i]
    sel = np.isin(det, Q.TRACKER) & (a["h_len"][i] > 0)
    if not sel.any():
        return []
    idx = np.nonzero(sel)[0]
    t0, x0, y0, z0 = a["h_t0"][i][idx], a["h_x0"][i][idx], a["h_y0"][i][idx], a["h_z0"][i][idx]
    x1, y1, z1 = a["h_x1"][i][idx], a["h_y1"][i][idx], a["h_z1"][i][idx]
    hq, vol, lay = a["h_q"][i][idx], a["h_vol"][i][idx], a["h_layer"][i][idx]
    dE, ln, dd = a["h_dE"][i][idx], a["h_len"][i][idx], det[idx]
    hits = {}
    for k in np.argsort(t0):  # first crossing of each sensor per quirk
        key = (int(vol[k]), int(hq[k]))
        if key in hits:
            continue
        xm, ym, zm = (x0[k] + x1[k]) / 2, (y0[k] + y1[k]) / 2, (z0[k] + z1[k]) / 2
        late = t0[k] - math.sqrt(xm * xm + ym * ym + zm * zm) / Q.C_MM_NS > Q.TOF_WINDOW
        hits[key] = dict(x=xm, y=ym, z=zm, dE=dE[k], len=ln[k], det=int(dd[k]), lay=int(lay[k]), vol=int(vol[k]),
                         q={int(hq[k])}, late=late)
    out = []
    for v in sorted({k[0] for k in hits}):
        h = [hits[(v, q)] for q in (0, 1) if (v, q) in hits]
        if len(h) == 2:
            a0, b0 = h
            p1, p2 = Q.PITCH[a0["det"]]
            r = math.hypot(a0["x"], a0["y"])
            drphi = r * Q.dphi(math.atan2(a0["y"], a0["x"]), math.atan2(b0["y"], b0["x"]))
            d2 = abs(a0["z"] - b0["z"]) if a0["det"] in (Q.PXB, Q.TIB, Q.TOB) else abs(math.hypot(a0["x"], a0["y"]) - math.hypot(b0["x"], b0["y"]))
            if drphi < 2 * p1 * MERGE and d2 < (2 * p2 * MERGE if a0["det"] in Q.PIXEL else p2):
                e = a0["dE"] + b0["dE"]
                m = dict(a0)
                m.update(x=(a0["x"] * a0["dE"] + b0["x"] * b0["dE"]) / e, y=(a0["y"] * a0["dE"] + b0["y"] * b0["dE"]) / e,
                         z=(a0["z"] + b0["z"]) / 2, dE=e, len=min(a0["len"], b0["len"]), q={0, 1}, late=a0["late"] and b0["late"])
                h = [m]
        for c in h:
            p1, _ = Q.PITCH[c["det"]]
            c["srphi"] = 0.010 if c["det"] in Q.PIXEL else p1 / math.sqrt(12.)
            c["szz"] = 0.020 if c["det"] == Q.PXB else (Q.ZRES_STEREO[c["det"]] if (c["det"], c["lay"]) in Q.STEREO else 0.)
            c["nmeas"] = 2 if (c["det"], c["lay"]) in Q.STEREO else 1
            th = int(thick[c["vol"]])
            x = c["dE"] / (c["len"] / 10.)
            x *= 1 + rng.normal(0., Q.NOISE["pix" if c["det"] in Q.PIXEL else "strip"])
            c["dedx"], c["p"] = x, mip.cdf(th, x)
            out.append(c)
    return out


def candidate(cl, qset, vz, beta, p, rng):
    """Track candidate from the clusters containing the quirks in qset."""
    mine = [c for c in cl if c["q"] & qset]
    valid = [c for c in mine if not c["late"]]
    if len(valid) < 3:
        return None
    ft = Q.fit_track(valid, vz, beta * p, rng, beta)
    if ft is None:
        return None
    pix = [c for c in valid if c["det"] in Q.PIXEL and not (c["det"] == Q.PXB and c["lay"] == 1)]
    strip = [c for c in valid if c["det"] not in Q.PIXEL]
    ft.update(npix=len(pix), fvalid=len(valid) / len(mine), ndedx=len(pix) + sum(c["nmeas"] for c in strip),
              gstrip=Q.gstrip(np.array([c["p"] for c in strip])), fpix=Q.fpix(np.array([c["p"] for c in pix])),
              ih=Q.harmonic2([c["dedx"] for c in pix + strip]), nmerged=sum(len(c["q"]) == 2 for c in valid), nhit=len(valid))
    return ft


def muon_obj(a, i, q):
    """First muon-station hit of quirk q: (n stations, beta_eff, |eta|, phi) or None."""
    det = a["h_det"][i]
    m = np.isin(det, (Q.MB, Q.ME)) & (a["h_q"][i] == q)
    if not m.any():
        return None
    st = {(int(d), int(l)) for d, l in zip(det[m], a["h_layer"][i][m])}
    k = np.argmin(np.where(m, a["h_t0"][i], 1e30))
    x, y, z, t = a["h_x0"][i][k], a["h_y0"][i][k], a["h_z0"][i][k], a["h_t0"][i][k]
    return len(st), math.sqrt(x * x + y * y + z * z) / (Q.C_MM_NS * t), abs(math.asinh(z / math.hypot(x, y))), math.atan2(y, x)


def trigger(a, i):
    """Event trigger efficiency (nominal, up, down): stau maps vs the effective beta to the muon stations."""
    objs = [o for o in (muon_obj(a, i, 0), muon_obj(a, i, 1)) if o and o[0] >= 2]
    if len(objs) == 2 and math.hypot(objs[0][2] - objs[1][2], Q.dphi(objs[0][3], objs[1][3])) < 0.1:
        objs = [max(objs, key=lambda o: o[1])]  # coincident pair: one muon object
    return [1. - np.prod([1. - Q.trig_eff(min(o[1], 0.999), o[2], var) for o in objs]) for var in (0, 1, 2)]


def run(outdir, pattern, fout, seed=5):
    rng = np.random.default_rng(seed)
    mip = Q.MipTemplates()
    files = sorted(glob.glob("%s/%s.root" % (outdir, pattern)))
    rec, ngen, nkept = [], 0, 0
    for f in files:
        base = f[:-5]
        meta, thick = Q.load_meta(base + "_meta.json")
        g = np.load(base + ".npz")
        cols = {c: k for k, c in enumerate(g["cols"])}
        gd = {int(r[0]): r for r in g["data"]}
        ngen += int(g["nin"])
        nkept += len(g["data"])
        for a in Q.read_events(f):
            for i in range(len(a["evt"])):
                row = gd.get(int(a["genEvt"][i]))
                if row is None:
                    continue
                rec.append(event(a, i, row, cols, thick, mip, rng))
    keys = sorted(rec[0].keys()) if rec else []
    np.savez(fout, ngen=ngen, nkept=nkept, nsim=len(rec), keys=np.array(keys),
             data=np.array([[r[k] for k in keys] for r in rec], dtype=float), cuts=np.array(CUTS))
    print("emulate: %s %d files, gen %d, kept %d, sim %d" % (pattern, len(files), ngen, nkept, len(rec)))


def event(a, i, row, cols, thick, mip, rng):
    mq = row[cols["m"]]
    bq = [row[cols["q%d_p" % k]] / math.hypot(row[cols["q%d_p" % k]], mq) for k in (1, 2)]
    ppr, mpr = row[cols["pr_p"]], row[cols["pr_m"]]
    bpr = ppr / math.hypot(ppr, mpr)
    vz = row[cols["vz"]]
    cl = clusters(a, i, thick, mip, rng, bq)
    cands = []
    pq = [row[cols["q%d_p" % k]] for k in (1, 2)]
    c0 = candidate(cl, {0}, vz, bq[0], pq[0], rng)
    c1 = candidate(cl, {1}, vz, bq[1], pq[1], rng)
    if c0 is not None and c0["nmerged"] >= 0.5 * c0["nhit"]:  # bound pair: refit with the pair momentum
        c0 = candidate(cl, {0}, vz, bpr, ppr, rng)
    for c, qs in ((c0, 0), (c1, 1)):
        if c is None:
            continue
        merged = c["nmerged"] >= 0.5 * c["nhit"]
        if merged and qs == 1 and c0 is not None:
            continue  # same merged track as c0
        d = "pr" if merged else "q%d" % (qs + 1)
        c["gdir"] = d
        c["beta"] = bpr if merged else bq[qs]
        r = Q.mini_cone(c["pt"])
        mini = Q.cone_sum(row, cols, d, "all", r)
        tk = Q.cone_sum(row, cols, d, "ch", 0.3)
        if not merged:  # partner quirk inside the cones
            o = "q%d" % (2 - qs)
            dr = math.hypot(row[cols[d + "_eta"]] - row[cols[o + "_eta"]], Q.dphi(row[cols[d + "_phi"]], row[cols[o + "_phi"]]))
            ptp = row[cols[o + "_pt"]]
            mini += ptp if dr < r else 0.
            tk += ptp if dr < 0.3 else 0.
        c["miniiso"] = mini / c["pt"]
        c["tkiso"] = tk
        nd = len(a["q_edepDet"][i]) // 2
        qs_set = (0, 1) if merged else (qs,)
        ecal = sum(a["q_edepDet"][i][q * nd + k] for q in qs_set for k in Q.CALO) / 1000.  # quirk's own deposit, GeV
        c["eop"] = (row[cols[d + "_e0.3"]] + ecal) / c["p"]
        c["ecal"] = ecal
        cands.append(c)
    trig = trigger(a, i)
    # cumulative cutflow: the event passes step k if one candidate passes steps 2..k
    tests = [lambda c: c["pt"] > 55, lambda c: abs(c["eta"]) < 1, lambda c: c["npix"] >= 2, lambda c: c["fvalid"] > 0.8,
             lambda c: c["ndedx"] >= 10, lambda c: True, lambda c: c["chi2ndof"] < 5, lambda c: c["dz"] < 1.,
             lambda c: c["dxy"] < 0.2, lambda c: c["miniiso"] < 0.02, lambda c: c["tkiso"] < 15, lambda c: c["eop"] < 0.3,
             lambda c: c["sinvpt"] < 8e-4, lambda c: 0.3 < c["fpix"] <= 1.]
    best, depth = None, -1
    for c in cands:
        k = 0
        while k < len(tests) and tests[k](c):
            k += 1
        if k > depth or (k == depth and c["gstrip"] > best["gstrip"]):
            depth, best = k, c
    depth = max(depth, 0)
    r = dict(depth=depth, trig=trig[0], trig_up=trig[1], trig_dn=trig[2], ncand=len(cands),
             maxSep=float(a["maxSep"][i]), cpu=float(a["cpu"][i]))
    for k in ("ecal", "eop", "pt", "p", "eta", "gstrip", "fpix", "ih", "beta", "chi2ndof", "dxy", "dz", "sinvpt", "miniiso", "tkiso",
              "nhit", "nmerged", "npix", "ndedx", "fvalid"):
        r["c_" + k] = best[k] if best else -1.
    r["c_merged"] = 1. if best and best["gdir"] == "pr" else 0.
    return r


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], sys.argv[3])
