"""Emulated EXO-18-002 reconstruction on quirksim output (CMS-like card): clusters, track fit, dE/dx
discriminants (Gstrip, Fpix, Ih), isolation from the GEN summary, muon trigger from the HEPData maps."""
import json
import math
import os
import numpy as np
import uproot
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
HEPDATA = os.path.join(HERE, "..", "hepdata", "original")
C_MM_NS = 299.792458
BFIELD = 3.8

PXB, PXF, TIB, TOB, TID, TEC, MB, ME = 2, 4, 5, 6, 7, 8, 16, 18
TRACKER = (PXB, PXF, TIB, TOB, TID, TEC)
PIXEL = (PXB, PXF)
# rphi pitch, second-coordinate pitch (pixel z / r) or strip length [mm]
PITCH = {PXB: (0.10, 0.15), PXF: (0.10, 0.15), TIB: (0.10, 120.), TOB: (0.15, 190.), TID: (0.12, 100.), TEC: (0.15, 100.)}
STEREO = {(TIB, 1), (TIB, 2), (TOB, 1), (TOB, 2)}  # double-sided: 2 measurements, z from stereo
ZRES_STEREO = {TIB: 0.23, TOB: 0.53}
NOISE = {"pix": float(os.environ.get("QP_NOISE_PIX", 0.40)), "strip": float(os.environ.get("QP_NOISE_STRIP", 0.15))}  # relative per-cluster smearing
TAG = "p%g_s%g" % (NOISE["pix"], NOISE["strip"])
TOF_WINDOW = 25.  # ns, late hits are lost
HITERR = float(os.environ.get("QP_HITERR", 1.))  # hit resolution scale in the fit (systematic)
BRANCHES = ["evt", "genEvt", "q_pdg", "q_px0", "q_py0", "q_pz0", "q_e0", "h_x0", "h_y0", "h_z0", "h_t0", "h_x1",
            "h_y1", "h_z1", "h_t1", "h_dE", "h_len", "h_q", "h_vol", "h_det", "h_layer", "maxSep", "cpu", "q_edepDet"]
CALO = (9, 10, 11, 12, 13, 15)  # EB ES EE HB HE HO


def load_meta(path):
    m = json.load(open(path))
    thick = []
    for v in m["vols"]:
        barrel = (v["zmax"] - v["zmin"]) > (v["rmax"] - v["rmin"])
        thick.append(round(1000 * ((v["rmax"] - v["rmin"]) if barrel else (v["zmax"] - v["zmin"]))))
    return m, np.array(thick)


def read_events(root_path, step=2000):
    with uproot.open(root_path) as f:
        for arrs in f["events"].iterate(BRANCHES, step_size=step, library="np"):
            yield arrs


# ---------------------------------------------------------------- dE/dx

class MipTemplates:
    """Per sensor thickness: sorted MIP dE/dx (MeV/cm, noise included)."""

    def __init__(self, path=os.path.join(HERE, "mip_templates_%s.npz" % TAG)):
        d = np.load(path)
        self.t = {int(k[1:]): d[k] for k in d.files if k.startswith("t")}
        self.mpv = {k: float(np.median(v)) for k, v in self.t.items()}

    def cdf(self, thick, x):
        a = self.t.get(thick)
        if a is None:  # nearest measured thickness
            a = self.t[min(self.t, key=lambda k: abs(k - thick))]
        return np.searchsorted(a, x) / len(a)


def gstrip(p):
    """Ias-like discriminant from per-hit MIP CDF values."""
    p = np.sort(np.clip(p, 0., 1.))
    n = len(p)
    if n == 0:
        return -1.
    i = np.arange(1, n + 1)
    return 3. / n * (1. / (12 * n) + np.sum(p * (p - (2 * i - 1) / (2. * n)) ** 2))


def fpix(p):
    """1 - combined probability (Fisher) that a MIP gives more charge on every hit."""
    q = np.clip(1. - np.asarray(p), 1e-12, 1.)
    if len(q) == 0:
        return -1.
    lp = -np.sum(np.log(q))
    s, t = 0., 1.
    for k in range(len(q)):
        s += t
        t *= lp / (k + 1)
    return 1. - math.exp(-lp) * s


def harmonic2(x):
    x = np.asarray(x)
    return float(np.mean(x ** -2) ** -0.5) if len(x) else -1.


# ---------------------------------------------------------------- track fit

def circle_fit(x, y, w):
    """Weighted algebraic circle fit (Kasa); returns xc, yc, R."""
    A = np.c_[x, y, np.ones_like(x)]
    b = -(x * x + y * y)
    sw = np.sqrt(w)
    sol, *_ = np.linalg.lstsq(A * sw[:, None], b * sw, rcond=None)
    xc, yc = -sol[0] / 2, -sol[1] / 2
    r2 = xc * xc + yc * yc - sol[2]
    return xc, yc, math.sqrt(max(r2, 1e-9))


def fit_track(cl, vz, betap, rng, beta=1.):
    """cl: list of cluster dicts. Returns dict with pt, chi2ndof, dxy, dz, sig_invpt, eta, phi."""
    x = np.array([c["x"] for c in cl])
    y = np.array([c["y"] for c in cl])
    z = np.array([c["z"] for c in cl])
    r = np.hypot(x, y)
    srp = np.array([c["srphi"] for c in cl]) * HITERR
    # multiple scattering floor, about 0.4 X0 of tracker
    th0 = 13.6e-3 * 0.6 / max(betap, 1e-3)  # betap = beta p [GeV]
    srp_eff = np.sqrt(srp ** 2 + (th0 * r / math.sqrt(3.)) ** 2)
    phi = np.arctan2(y, x)
    xs = x + rng.normal(0., srp, len(x)) * -np.sin(phi)
    ys = y + rng.normal(0., srp, len(x)) * np.cos(phi)
    w = 1. / srp_eff ** 2
    if len(x) < 3:
        return None
    # include the beam line weakly so 3-point fits stay defined
    xc, yc, R = circle_fit(np.r_[xs, 0.], np.r_[ys, 0.], np.r_[w, 1e-2])
    res = np.hypot(xs - xc, ys - yc) - R
    chi2 = np.sum(res ** 2 * w)
    ndof = len(x) - 3
    dxy = abs(math.hypot(xc, yc) - R)
    # z line vs transverse radius, on pixel and stereo clusters
    zm = np.array([c["szz"] > 0 for c in cl])
    dz, tanl = 0., 0.
    if zm.sum() >= 2:
        sz = np.array([c["szz"] for c in cl])[zm]
        zz = z[zm] + rng.normal(0., sz)
        ww = 1. / (sz ** 2 + (th0 * r[zm] / math.sqrt(3.)) ** 2)
        A = np.c_[np.ones(zm.sum()), r[zm]]
        cov = np.linalg.inv(A.T @ (A * ww[:, None]))
        z0, tanl = cov @ (A.T @ (ww * zz))
        chi2 += np.sum((zz - z0 - tanl * r[zm]) ** 2 * ww)
        ndof += zm.sum() - 2
        dz = abs(z0 - vz)
    pt = 0.3 * BFIELD * R / 1000.  # GeV, R in mm
    # Gluckstern + multiple scattering for sigma(1/pT)
    L = (r.max() - r.min()) / 1000.
    n = len(x)
    s_rphi = float(np.mean(srp)) / 1000.
    sk = s_rphi / max(L, 0.05) ** 2 * math.sqrt(720. / (n + 4))
    sms = 0.0136 / (max(L, 0.05) * max(beta, 1e-3)) * math.sqrt(0.4) * 1.43 / (0.3 * BFIELD)
    sinvpt = math.hypot(sk / (0.3 * BFIELD), sms / max(pt, 1e-3) * 1.)
    # direction at the vertex: tangent of the circle at the point of closest approach
    phi0 = math.atan2(yc, xc) + (math.pi / 2 if (xs[0] - xc) * ys[0] - (ys[0] - yc) * xs[0] > 0 else -math.pi / 2)
    eta = math.asinh(tanl) if zm.sum() >= 2 else math.asinh(z[np.argmax(r)] / r.max())
    return dict(pt=pt, chi2ndof=chi2 / max(ndof, 1), dxy=dxy, dz=dz, sinvpt=sinvpt, eta=eta, phi=phi0, p=pt * math.cosh(eta))


# ---------------------------------------------------------------- trigger

def load_trigger():
    """Stau trigger efficiency vs beta in |eta| bins (App. Fig. 12): nominal, up, down."""
    names = ["top_left", "top_middle", "top_right", "bottom_left", "bottom_middle", "bottom_right"]
    edges = [0., 0.3, 0.6, 0.9, 1.2, 2.1, 2.4]
    tab = []
    for n in names:
        d = yaml.safe_load(open(os.path.join(HEPDATA, "appendix_figure_12_%s.yaml" % n)))
        tab.append([[float(v["value"]) for v in dv["values"]] for dv in d["dependent_variables"]])
    return np.array(edges), np.array(tab)  # tab[eta][nom/up/down][beta bin]


TRIG_EDGES, TRIG_TAB = None, None


def trig_eff(beta, aeta, var=0):
    global TRIG_EDGES, TRIG_TAB
    if TRIG_TAB is None:
        TRIG_EDGES, TRIG_TAB = load_trigger()
    if aeta >= 2.4 or beta <= 0:
        return 0.
    ie = np.searchsorted(TRIG_EDGES, aeta, side="right") - 1
    ib = min(int(beta * 10), 9)
    return TRIG_TAB[ie][var][ib]


# ---------------------------------------------------------------- isolation

CONES = (0.05, 0.1, 0.2, 0.3)


def cone_sum(row, cols, d, kind, r):
    v = [row[cols["%s_%s%g" % (d, kind, c)]] for c in CONES]
    return float(np.interp(r, CONES, v))


def mini_cone(pt):
    return 0.2 if pt < 50 else (10. / pt if pt < 200 else 0.05)


def dphi(a, b):
    d = abs(a - b) % (2 * math.pi)
    return 2 * math.pi - d if d > math.pi else d
