#!/usr/bin/env python3
"""GEN filter + summary. Keep HepMC3 events where a quirk or the pair has |eta| < etamax and write
per-event generator quantities (quirks, pair, isolation sums around q1, q2, pair).
usage: prep_gen.py in.hepmc3 out.hepmc3 out.npz [etamax=1.3] [bound]
bound: keep only pairs with |eta| < 1.1 and beta > 0.4 (bound-pair selection needs both)"""
import math
import sys
import numpy as np

fin, fout, fnpz = sys.argv[1:4]
etamax = float(sys.argv[4]) if len(sys.argv) > 4 else 1.3
bound = len(sys.argv) > 5 and sys.argv[5] == "bound"
CHG = {11, 13, 211, 321, 2212, 3112, 3222, 3312, 3334}  # stable charged after Pythia decays
INVIS = {12, 14, 16}
CONES = (0.05, 0.1, 0.2, 0.3)


def kin(px, py, pz):
    pt = math.hypot(px, py)
    p = math.sqrt(pt * pt + pz * pz)
    eta = math.asinh(pz / pt) if pt > 0 else math.copysign(99., pz)
    return pt, eta, math.atan2(py, px), p


def dr2(e1, p1, e2, p2):
    dp = abs(p1 - p2)
    dp = 2 * math.pi - dp if dp > math.pi else dp
    return (e1 - e2) ** 2 + dp * dp


rows, nin = [], 0


def process(ev, o):
    global nin
    nin += 1
    q, vis, vtx = [], [], {}
    for l in ev:
        if l.startswith("V "):
            f = l.split()
            if "@" in f:
                vtx[int(f[1])] = float(f[f.index("@") + 3])
        elif l.startswith("P ") and l.split()[-1] == "1":
            f = l.split()
            pdg = int(f[3])
            px, py, pz, e = map(float, f[4:8])
            if abs(pdg) == 17:
                q.append((pdg, px, py, pz, e, float(f[8]), int(f[2])))
            elif abs(pdg) not in INVIS:
                pt, eta, phi, _ = kin(px, py, pz)
                vis.append((pt, eta, phi, abs(pdg) in CHG, e))
    if len(q) != 2:
        return
    q.sort(key=lambda x: -x[0])  # q1 = positive PDG
    k1 = kin(*q[0][1:4])
    k2 = kin(*q[1][1:4])
    kp = kin(q[0][1] + q[1][1], q[0][2] + q[1][2], q[0][3] + q[1][3])
    if min(abs(k1[1]), abs(k2[1]), abs(kp[1])) > etamax:
        return
    if bound:
        ep = q[0][4] + q[1][4]
        if abs(kp[1]) > 1.1 or kp[3] / ep < 0.4:
            return
    iso = []
    for d in (k1, k2, kp):
        ch = [0.] * len(CONES)
        al = [0.] * len(CONES)
        en = 0.
        for pt, eta, phi, c, e in vis:
            r2 = dr2(d[1], d[2], eta, phi)
            for i, r in enumerate(CONES):
                if r2 < r * r:
                    al[i] += pt
                    if c and pt > 1.:
                        ch[i] += pt
            if r2 < 0.09:
                en += e
        iso.append(ch + al + [en])
    mq = q[0][5]
    mpair = math.sqrt(max((q[0][4] + q[1][4]) ** 2 - sum((q[0][i] + q[1][i]) ** 2 for i in (1, 2, 3)), 0.))
    genEvt = int(ev[0].split()[1])
    rows.append([genEvt, mq, *k1, *k2, *kp, mpair, vtx.get(q[0][6], 0.)] + sum(iso, []))
    for l in ev:
        o.write(l)


with open(fin) as f, open(fout, "w") as o:
    ev = None
    for l in f:
        if l.startswith("E "):
            if ev:
                process(ev, o)
            ev = [l]
        elif l.startswith("HepMC::Asciiv3-END"):
            if ev:
                process(ev, o)
            ev = None
        elif ev is None:
            o.write(l)
        else:
            ev.append(l)
    o.write("HepMC::Asciiv3-END_EVENT_LISTING\n")

cols = ["genEvt", "m", "q1_pt", "q1_eta", "q1_phi", "q1_p", "q2_pt", "q2_eta", "q2_phi", "q2_p",
        "pr_pt", "pr_eta", "pr_phi", "pr_p", "pr_m", "vz"]
for d in ("q1", "q2", "pr"):
    cols += ["%s_ch%g" % (d, r) for r in CONES] + ["%s_all%g" % (d, r) for r in CONES] + ["%s_e0.3" % d]
np.savez(fnpz, data=np.array(rows, dtype=float), cols=np.array(cols), nin=nin, etamax=etamax)
print("prep_gen: kept %d of %d" % (len(rows), nin))
