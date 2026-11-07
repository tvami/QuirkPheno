#!/usr/bin/env python3
"""EXO-18-002 selection emulated on quirksim output. usage: emulate.py <outdir> <tag-glob> <out.npz>
Writes the per-event cutflow and the best candidate (pt, p, Gstrip, Fpix, Ih) for every trigger variation."""
import glob
import math
import sys
import numpy as np
import os
import qreco

MERGE = float(os.environ.get("QP_MERGE", 1.))  # two-hit resolution scale (systematic)

CUTS = ["All", "Trigger", "pT>55", "|eta|<1", "Npix>=2", "fvalid>0.8", "NdEdx>=10", "HighPurity", "chi2<5",
        "dz<0.1", "dxy<0.02", "MiniIso<0.02", "TkIso<15", "E/p<0.3", "sigpt<8e-4", "Fpix>0.3"]


def clusters(arrays, ievt, thick, mip, rng, beta_q):
    """Tracker clusters of one event: per (vol) merge the two quirks' hits within the two-hit resolution."""
    det = arrays["h_det"][ievt]
    sel = np.isin(det, qreco.TRACKER) & (arrays["h_len"][ievt] > 0)
    if not sel.any():
        return []
    idx = np.nonzero(sel)[0]
    t0, x0, y0, z0 = arrays["h_t0"][ievt][idx], arrays["h_x0"][ievt][idx], arrays["h_y0"][ievt][idx], arrays["h_z0"][ievt][idx]
    x1, y1, z1 = arrays["h_x1"][ievt][idx], arrays["h_y1"][ievt][idx], arrays["h_z1"][ievt][idx]
    hq, vol, lay = arrays["h_q"][ievt][idx], arrays["h_vol"][ievt][idx], arrays["h_layer"][ievt][idx]
    dE, ln, dd = arrays["h_dE"][ievt][idx], arrays["h_len"][ievt][idx], det[idx]
    hits = {}
    for ihit in np.argsort(t0):  # first crossing of each sensor per quirk
        key = (int(vol[ihit]), int(hq[ihit]))
        if key in hits:
            continue
        xm, ym, zm = (x0[ihit] + x1[ihit]) / 2, (y0[ihit] + y1[ihit]) / 2, (z0[ihit] + z1[ihit]) / 2
        late = t0[ihit] - math.sqrt(xm * xm + ym * ym + zm * zm) / qreco.C_MM_NS > qreco.TOF_WINDOW
        hits[key] = dict(x=xm, y=ym, z=zm, dE=dE[ihit], len=ln[ihit], det=int(dd[ihit]), lay=int(lay[ihit]), vol=int(vol[ihit]),
                         q={int(hq[ihit])}, late=late)
    out = []
    for volume in sorted({key[0] for key in hits}):
        group = [hits[(volume, quirk)] for quirk in (0, 1) if (volume, quirk) in hits]
        if len(group) == 2:
            a0, b0 = group
            p1, p2 = qreco.PITCH[a0["det"]]
            rad = math.hypot(a0["x"], a0["y"])
            drphi = rad * qreco.dphi(math.atan2(a0["y"], a0["x"]), math.atan2(b0["y"], b0["x"]))
            d2 = abs(a0["z"] - b0["z"]) if a0["det"] in (qreco.PXB, qreco.TIB, qreco.TOB) \
                else abs(math.hypot(a0["x"], a0["y"]) - math.hypot(b0["x"], b0["y"]))
            if drphi < 2 * p1 * MERGE and d2 < (2 * p2 * MERGE if a0["det"] in qreco.PIXEL else p2):
                esum = a0["dE"] + b0["dE"]
                merged = dict(a0)
                merged.update(x=(a0["x"] * a0["dE"] + b0["x"] * b0["dE"]) / esum,
                              y=(a0["y"] * a0["dE"] + b0["y"] * b0["dE"]) / esum,
                              z=(a0["z"] + b0["z"]) / 2, dE=esum, len=min(a0["len"], b0["len"]), q={0, 1},
                              late=a0["late"] and b0["late"])
                group = [merged]
        for clu in group:
            p1, _ = qreco.PITCH[clu["det"]]
            clu["srphi"] = 0.010 if clu["det"] in qreco.PIXEL else p1 / math.sqrt(12.)
            clu["szz"] = 0.020 if clu["det"] == qreco.PXB else (
                qreco.ZRES_STEREO[clu["det"]] if (clu["det"], clu["lay"]) in qreco.STEREO else 0.)
            clu["nmeas"] = 2 if (clu["det"], clu["lay"]) in qreco.STEREO else 1
            th = int(thick[clu["vol"]])
            dedx = clu["dE"] / (clu["len"] / 10.)
            dedx *= 1 + rng.normal(0., qreco.NOISE["pix" if clu["det"] in qreco.PIXEL else "strip"])
            clu["dedx"], clu["p"] = dedx, mip.cdf(th, dedx)
            out.append(clu)
    return out


def candidate(cl, qset, vz, beta, mom, rng):
    """Track candidate from the clusters containing the quirks in qset."""
    mine = [clu for clu in cl if clu["q"] & qset]
    valid = [clu for clu in mine if not clu["late"]]
    if len(valid) < 3:
        return None
    ft = qreco.fit_track(valid, vz, beta * mom, rng, beta)
    if ft is None:
        return None
    pix = [clu for clu in valid if clu["det"] in qreco.PIXEL and not (clu["det"] == qreco.PXB and clu["lay"] == 1)]
    strip = [clu for clu in valid if clu["det"] not in qreco.PIXEL]
    ft.update(npix=len(pix), fvalid=len(valid) / len(mine), ndedx=len(pix) + sum(clu["nmeas"] for clu in strip),
              gstrip=qreco.gstrip(np.array([clu["p"] for clu in strip])), fpix=qreco.fpix(np.array([clu["p"] for clu in pix])),
              ih=qreco.harmonic2([clu["dedx"] for clu in pix + strip]), nmerged=sum(len(clu["q"]) == 2 for clu in valid),
              nhit=len(valid))
    return ft


def muon_obj(arrays, ievt, quirk):
    """First muon-station hit of quirk q: (n stations, beta_eff, |eta|, phi) or None."""
    det = arrays["h_det"][ievt]
    mask = np.isin(det, (qreco.MB, qreco.ME)) & (arrays["h_q"][ievt] == quirk)
    if not mask.any():
        return None
    st = {(int(det_id), int(layer)) for det_id, layer in zip(det[mask], arrays["h_layer"][ievt][mask])}
    first = np.argmin(np.where(mask, arrays["h_t0"][ievt], 1e30))
    hx, hy, hz, ht = (arrays["h_x0"][ievt][first], arrays["h_y0"][ievt][first], arrays["h_z0"][ievt][first],
                      arrays["h_t0"][ievt][first])
    return (len(st), math.sqrt(hx * hx + hy * hy + hz * hz) / (qreco.C_MM_NS * ht),
            abs(math.asinh(hz / math.hypot(hx, hy))), math.atan2(hy, hx))


def trigger(arrays, ievt):
    """Event trigger efficiency (nominal, up, down): stau maps vs the effective beta to the muon stations."""
    objs = [obj for obj in (muon_obj(arrays, ievt, 0), muon_obj(arrays, ievt, 1)) if obj and obj[0] >= 2]
    if len(objs) == 2 and math.hypot(objs[0][2] - objs[1][2], qreco.dphi(objs[0][3], objs[1][3])) < 0.1:
        objs = [max(objs, key=lambda obj: obj[1])]  # coincident pair: one muon object
    return [1. - np.prod([1. - qreco.trig_eff(min(obj[1], 0.999), obj[2], var) for obj in objs]) for var in (0, 1, 2)]


def run(outdir, pattern, fout, seed=5):
    rng = np.random.default_rng(seed)
    mip = qreco.MipTemplates()
    files = sorted(glob.glob("%s/%s.root" % (outdir, pattern)))
    rec, ngen, nkept = [], 0, 0
    for fname in files:
        base = fname[:-5]
        meta, thick = qreco.load_meta(base + "_meta.json")
        gen = np.load(base + ".npz")
        cols = {name: col for col, name in enumerate(gen["cols"])}
        gd = {int(grow[0]): grow for grow in gen["data"]}
        ngen += int(gen["nin"])
        nkept += len(gen["data"])
        for arrays in qreco.read_events(fname):
            for ievt in range(len(arrays["evt"])):
                row = gd.get(int(arrays["genEvt"][ievt]))
                if row is None:
                    continue
                rec.append(event(arrays, ievt, row, cols, thick, mip, rng))
    keys = sorted(rec[0].keys()) if rec else []
    np.savez(fout, ngen=ngen, nkept=nkept, nsim=len(rec), keys=np.array(keys),
             data=np.array([[res[key] for key in keys] for res in rec], dtype=float), cuts=np.array(CUTS))
    print("emulate: %s %d files, gen %d, kept %d, sim %d" % (pattern, len(files), ngen, nkept, len(rec)))


def event(arrays, ievt, row, cols, thick, mip, rng):
    mq = row[cols["m"]]
    bq = [row[cols["q%d_p" % qn]] / math.hypot(row[cols["q%d_p" % qn]], mq) for qn in (1, 2)]
    ppr, mpr = row[cols["pr_p"]], row[cols["pr_m"]]
    bpr = ppr / math.hypot(ppr, mpr)
    vz = row[cols["vz"]]
    cl = clusters(arrays, ievt, thick, mip, rng, bq)
    cands = []
    pq = [row[cols["q%d_p" % qn]] for qn in (1, 2)]
    c0 = candidate(cl, {0}, vz, bq[0], pq[0], rng)
    c1 = candidate(cl, {1}, vz, bq[1], pq[1], rng)
    if c0 is not None and c0["nmerged"] >= 0.5 * c0["nhit"]:  # bound pair: refit with the pair momentum
        c0 = candidate(cl, {0}, vz, bpr, ppr, rng)
    for cand, qs in ((c0, 0), (c1, 1)):
        if cand is None:
            continue
        merged = cand["nmerged"] >= 0.5 * cand["nhit"]
        if merged and qs == 1 and c0 is not None:
            continue  # same merged track as c0
        gdir = "pr" if merged else "q%d" % (qs + 1)
        cand["gdir"] = gdir
        cand["beta"] = bpr if merged else bq[qs]
        cone = qreco.mini_cone(cand["pt"])
        mini = qreco.cone_sum(row, cols, gdir, "all", cone)
        tk = qreco.cone_sum(row, cols, gdir, "ch", 0.3)
        if not merged:  # partner quirk inside the cones
            other = "q%d" % (2 - qs)
            dr = math.hypot(row[cols[gdir + "_eta"]] - row[cols[other + "_eta"]],
                            qreco.dphi(row[cols[gdir + "_phi"]], row[cols[other + "_phi"]]))
            ptp = row[cols[other + "_pt"]]
            mini += ptp if dr < cone else 0.
            tk += ptp if dr < 0.3 else 0.
        cand["miniiso"] = mini / cand["pt"]
        cand["tkiso"] = tk
        nd = len(arrays["q_edepDet"][ievt]) // 2
        qs_set = (0, 1) if merged else (qs,)
        ecal = sum(arrays["q_edepDet"][ievt][quirk * nd + calo]
                   for quirk in qs_set for calo in qreco.CALO) / 1000.  # quirk's own deposit, GeV
        cand["eop"] = (row[cols[gdir + "_e0.3"]] + ecal) / cand["p"]
        cand["ecal"] = ecal
        cands.append(cand)
    trig = trigger(arrays, ievt)
    # cumulative cutflow: the event passes step k if one candidate passes steps 2..k
    tests = [lambda cand: cand["pt"] > 55, lambda cand: abs(cand["eta"]) < 1, lambda cand: cand["npix"] >= 2,
             lambda cand: cand["fvalid"] > 0.8, lambda cand: cand["ndedx"] >= 10, lambda cand: True,
             lambda cand: cand["chi2ndof"] < 5, lambda cand: cand["dz"] < 1., lambda cand: cand["dxy"] < 0.2,
             lambda cand: cand["miniiso"] < 0.02, lambda cand: cand["tkiso"] < 15, lambda cand: cand["eop"] < 0.3,
             lambda cand: cand["sinvpt"] < 8e-4, lambda cand: 0.3 < cand["fpix"] <= 1.]
    best, depth = None, -1
    for cand in cands:
        step = 0
        while step < len(tests) and tests[step](cand):
            step += 1
        if step > depth or (step == depth and cand["gstrip"] > best["gstrip"]):
            depth, best = step, cand
    depth = max(depth, 0)
    result = dict(depth=depth, trig=trig[0], trig_up=trig[1], trig_dn=trig[2], ncand=len(cands),
                  maxSep=float(arrays["maxSep"][ievt]), cpu=float(arrays["cpu"][ievt]))
    for key in ("ecal", "eop", "pt", "p", "eta", "gstrip", "fpix", "ih", "beta", "chi2ndof", "dxy", "dz", "sinvpt", "miniiso",
                "tkiso", "nhit", "nmerged", "npix", "ndedx", "fvalid"):
        result["c_" + key] = best[key] if best else -1.
    result["c_merged"] = 1. if best and best["gdir"] == "pr" else 0.
    return result


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], sys.argv[3])
