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
    meta = json.load(open(path))
    thick = []
    for vol in meta["vols"]:
        barrel = (vol["zmax"] - vol["zmin"]) > (vol["rmax"] - vol["rmin"])
        thick.append(round(1000 * ((vol["rmax"] - vol["rmin"]) if barrel else (vol["zmax"] - vol["zmin"]))))
    return meta, np.array(thick)


def read_events(root_path, step=2000):
    with uproot.open(root_path) as handle:
        for arrs in handle["events"].iterate(BRANCHES, step_size=step, library="np"):
            yield arrs


# ---------------------------------------------------------------- dE/dx

class MipTemplates:
    """Per sensor thickness: sorted MIP dE/dx (MeV/cm, noise included)."""

    def __init__(self, path=os.path.join(HERE, "mip_templates_%s.npz" % TAG)):
        npz = np.load(path)
        self.t = {int(name[1:]): npz[name] for name in npz.files if name.startswith("t")}
        self.mpv = {thk: float(np.median(vals)) for thk, vals in self.t.items()}

    def cdf(self, thick, value):
        table = self.t.get(thick)
        if table is None:  # nearest measured thickness
            table = self.t[min(self.t, key=lambda thk: abs(thk - thick))]
        return np.searchsorted(table, value) / len(table)


def gstrip(probs):
    """Ias-like discriminant from per-hit MIP CDF values."""
    probs = np.sort(np.clip(probs, 0., 1.))
    nhit = len(probs)
    if nhit == 0:
        return -1.
    rank = np.arange(1, nhit + 1)
    return 3. / nhit * (1. / (12 * nhit) + np.sum(probs * (probs - (2 * rank - 1) / (2. * nhit)) ** 2))


def fpix(probs):
    """1 - combined probability (Fisher) that a MIP gives more charge on every hit."""
    comp = np.clip(1. - np.asarray(probs), 1e-12, 1.)
    if len(comp) == 0:
        return -1.
    lp = -np.sum(np.log(comp))
    total, term = 0., 1.
    for idx in range(len(comp)):
        total += term
        term *= lp / (idx + 1)
    return 1. - math.exp(-lp) * total


def harmonic2(vals):
    vals = np.asarray(vals)
    return float(np.mean(vals ** -2) ** -0.5) if len(vals) else -1.


# ---------------------------------------------------------------- track fit

def circle_fit(xpos, ypos, weight):
    """Weighted algebraic circle fit (Kasa); returns xc, yc, R."""
    design = np.c_[xpos, ypos, np.ones_like(xpos)]
    rhs = -(xpos * xpos + ypos * ypos)
    sw = np.sqrt(weight)
    sol, *_ = np.linalg.lstsq(design * sw[:, None], rhs * sw, rcond=None)
    xc, yc = -sol[0] / 2, -sol[1] / 2
    r2 = xc * xc + yc * yc - sol[2]
    return xc, yc, math.sqrt(max(r2, 1e-9))


def fit_track(cl, vz, betap, rng, beta=1.):
    """cl: list of cluster dicts. Returns dict with pt, chi2ndof, dxy, dz, sig_invpt, eta, phi."""
    xpos = np.array([clu["x"] for clu in cl])
    ypos = np.array([clu["y"] for clu in cl])
    zpos = np.array([clu["z"] for clu in cl])
    rad = np.hypot(xpos, ypos)
    srp = np.array([clu["srphi"] for clu in cl]) * HITERR
    # multiple scattering floor, about 0.4 X0 of tracker
    th0 = 13.6e-3 * 0.6 / max(betap, 1e-3)  # betap = beta p [GeV]
    srp_eff = np.sqrt(srp ** 2 + (th0 * rad / math.sqrt(3.)) ** 2)
    phi = np.arctan2(ypos, xpos)
    xs = xpos + rng.normal(0., srp, len(xpos)) * -np.sin(phi)
    ys = ypos + rng.normal(0., srp, len(xpos)) * np.cos(phi)
    weight = 1. / srp_eff ** 2
    if len(xpos) < 3:
        return None
    # include the beam line weakly so 3-point fits stay defined
    xc, yc, radius = circle_fit(np.r_[xs, 0.], np.r_[ys, 0.], np.r_[weight, 1e-2])
    res = np.hypot(xs - xc, ys - yc) - radius
    chi2 = np.sum(res ** 2 * weight)
    ndof = len(xpos) - 3
    dxy = abs(math.hypot(xc, yc) - radius)
    # z line vs transverse radius, on pixel and stereo clusters
    zm = np.array([clu["szz"] > 0 for clu in cl])
    dz, tanl = 0., 0.
    if zm.sum() >= 2:
        sz = np.array([clu["szz"] for clu in cl])[zm]
        zz = zpos[zm] + rng.normal(0., sz)
        ww = 1. / (sz ** 2 + (th0 * rad[zm] / math.sqrt(3.)) ** 2)
        design = np.c_[np.ones(zm.sum()), rad[zm]]
        cov = np.linalg.inv(design.T @ (design * ww[:, None]))
        z0, tanl = cov @ (design.T @ (ww * zz))
        chi2 += np.sum((zz - z0 - tanl * rad[zm]) ** 2 * ww)
        ndof += zm.sum() - 2
        dz = abs(z0 - vz)
    pt = 0.3 * BFIELD * radius / 1000.  # GeV, radius in mm
    # Gluckstern + multiple scattering for sigma(1/pT)
    lever = (rad.max() - rad.min()) / 1000.
    nhit = len(xpos)
    s_rphi = float(np.mean(srp)) / 1000.
    sk = s_rphi / max(lever, 0.05) ** 2 * math.sqrt(720. / (nhit + 4))
    sms = 0.0136 / (max(lever, 0.05) * max(beta, 1e-3)) * math.sqrt(0.4) * 1.43 / (0.3 * BFIELD)
    sinvpt = math.hypot(sk / (0.3 * BFIELD), sms / max(pt, 1e-3) * 1.)
    # direction at the vertex: tangent of the circle at the point of closest approach
    phi0 = math.atan2(yc, xc) + (math.pi / 2 if (xs[0] - xc) * ys[0] - (ys[0] - yc) * xs[0] > 0 else -math.pi / 2)
    eta = math.asinh(tanl) if zm.sum() >= 2 else math.asinh(zpos[np.argmax(rad)] / rad.max())
    return dict(pt=pt, chi2ndof=chi2 / max(ndof, 1), dxy=dxy, dz=dz, sinvpt=sinvpt, eta=eta, phi=phi0, p=pt * math.cosh(eta))


# ---------------------------------------------------------------- trigger

def load_trigger():
    """Stau trigger efficiency vs beta in |eta| bins (App. Fig. 12): nominal, up, down."""
    names = ["top_left", "top_middle", "top_right", "bottom_left", "bottom_middle", "bottom_right"]
    edges = [0., 0.3, 0.6, 0.9, 1.2, 2.1, 2.4]
    tab = []
    for name in names:
        doc = yaml.safe_load(open(os.path.join(HEPDATA, "appendix_figure_12_%s.yaml" % name)))
        tab.append([[float(val["value"]) for val in dv["values"]] for dv in doc["dependent_variables"]])
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


def cone_sum(row, cols, gdir, kind, radius):
    vals = [row[cols["%s_%s%g" % (gdir, kind, cone)]] for cone in CONES]
    return float(np.interp(radius, CONES, vals))


def mini_cone(pt):
    return 0.2 if pt < 50 else (10. / pt if pt < 200 else 0.05)


def dphi(phi1, phi2):
    diff = abs(phi1 - phi2) % (2 * math.pi)
    return 2 * math.pi - diff if diff > math.pi else diff
