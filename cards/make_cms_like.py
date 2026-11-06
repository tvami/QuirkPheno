#!/usr/bin/env python3
"""Writes the generic CMS-like detector card (cylinders and disks, all numbers public).

Sources (rounded; positions of the outer systems are approximate, see cards/README.md):
  [D] CMS Collaboration, "The CMS experiment at the CERN LHC", JINST 3 (2008) S08004
  [P] CMS Collaboration, "CMS Technical Design Report for the Pixel Detector Upgrade", CERN-LHCC-2012-016 (CMS-TDR-011)
  [M] tracker material budget vs eta, [D] Fig. 3.4 (about 0.4 X0 at eta = 0, peak about 1.8 X0 near |eta| = 1.4)

Options give the variations used to rank the simplifications (REPORT.md, question 3).
"""
import argparse

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--out", default="cms_like.card")
ap.add_argument("--tk-scale", type=float, default=1.0, help="scale of the passive tracker material (X/X0)")
ap.add_argument("--calo-scale", type=float, default=1.0, help="scale of the calorimeter absorber thickness")
ap.add_argument("--no-endcaps", action="store_true", help="barrel only")
ap.add_argument("--field", choices=["yoke", "uniform", "none"], default="yoke",
                help="yoke: 3.8 T in the solenoid and -1.8 T in the barrel iron; uniform: solenoid only")
a = ap.parse_args()

L = []
w = L.append


def tube(name, mat, rmin, rmax, zmin, zmax, det, layer=0, sens=False, calo=False, mirror=False):
    opts = "det=%s" % det + (" layer=%d" % layer if layer else "") + (" sens" if sens else "") + \
           (" calo" if calo else "") + (" mirror" if mirror else "")
    f = lambda v: v if isinstance(v, str) else "%.4f" % v
    w("tube %-10s %-28s %9s %11s %9s %11s %s" % (name, mat, f(rmin), f(rmax), f(zmin), f(zmax), opts))


tag = "tk%.2f_calo%.2f_%s_%s" % (a.tk_scale, a.calo_scale, "barrel" if a.no_endcaps else "full", a.field)
w("# Generic CMS-like detector, written by make_cms_like.py (%s)" % tag)
w("# lengths in cm; '+x0:v' / '+lam:v' / '+cm:v' / '+um:v' = thickness from the inner edge")
w("world 900 1150 G4_AIR")
w("checkOverlaps 1")
w("mix TkSupport 1.8 G4_C 0.55 G4_Al 0.15 G4_Cu 0.20 G4_KAPTON 0.10  # carbon fiber, cooling, cables")
w("mix TkServ 3.6 G4_C 0.55 G4_Al 0.15 G4_Cu 0.20 G4_KAPTON 0.10      # same, denser: cables and pipes in the TIB-TID gap")
w("mix MuonGas 0.00184 G4_Ar 0.837 G4_CARBON_DIOXIDE 0.163           # Ar/CO2 85/15 [D]")
w("mix CoilMix 2.81 G4_Al 0.86 G4_Cu 0.08 G4_Ti 0.03 G4_Nb 0.03       # Al-stabilized NbTi conductor and mandrel")

# ---------------- beam pipe [D,P]: Be, r = 2.2 cm, 0.8 mm
w("# beam pipe")
tube("PIPE_VAC", "G4_Galactic", 0, 2.17, -1100, 1100, "PIPE")
tube("PIPE", "G4_Be", 2.17, "+cm:0.08", -300, 300, "PIPE")
tube("PIPEFWD", "G4_STAINLESS-STEEL", 2.17, "+cm:0.15", 300.1, 1100, "PIPE", mirror=True)  # approximate

# ---------------- pixel [P]: 4 barrel layers, 3 disks per side, 285 um sensors
w("# pixel barrel and forward disks [P]; support per layer tuned to about 1.2% X0")
s = a.tk_scale
for i, r in enumerate([2.9, 6.8, 10.9, 16.0], 1):
    tube("PXB%d" % i, "G4_Si", r, "+um:285", -26.7, 26.7, "PXB", i, sens=True)
    tube("PXB%dS" % i, "TkSupport", r + 0.05, "+x0:%.5f" % (0.009 * s), -26.7, 26.7, "TKSUP", i)
if not a.no_endcaps:
    for i, z in enumerate([29.1, 39.6, 51.6], 1):
        tube("PXF%d" % i, "G4_Si", 4.5, 16.1, z, "+um:285", "PXF", i, sens=True, mirror=True)
        tube("PXF%dS" % i, "TkSupport", 4.5, 16.1, z + 0.05, "+x0:%.5f" % (0.012 * s), "TKSUP", i, mirror=True)
    tube("PIXSERV", "TkSupport", 4.5, 17.0, 56.0, "+x0:%.5f" % (0.10 * s), "TKSUP", mirror=True)

# ---------------- strips [D]: TIB 4, TOB 6 layers; TID 3, TEC 9 disks per side
w("# strip tracker [D]; layers 1-2 of TIB and TOB are double sided (two sensors)")
tib = [25.5, 33.9, 41.85, 49.8]
tob = [60.8, 69.2, 78.0, 86.8, 96.5, 108.0]
for i, r in enumerate(tib, 1):
    tube("TIB%d" % i, "G4_Si", r, "+um:%d" % (640 if i <= 2 else 320), -70, 70, "TIB", i, sens=True)
    tube("TIB%dS" % i, "TkSupport", r + 0.1, "+x0:%.5f" % (0.030 * s), -70, 70, "TKSUP", i)
for i, r in enumerate(tob, 1):
    tube("TOB%d" % i, "G4_Si", r, "+um:%d" % (1000 if i <= 2 else 500), -109, 109, "TOB", i, sens=True)
    tube("TOB%dS" % i, "TkSupport", r + 0.15, "+x0:%.5f" % (0.020 * s), -109, 109, "TKSUP", i)
tube("TKTUBE", "TkSupport", 116.0, "+x0:%.5f" % (0.020 * s), -280, 280, "TKSUP")
if not a.no_endcaps:
    for i, z in enumerate([80.0, 85.0, 90.0], 1):
        tube("TID%d" % i, "G4_Si", 20.0, 50.0, z, "+um:320", "TID", i, sens=True, mirror=True)
        tube("TID%dS" % i, "TkSupport", 20.0, 50.0, z + 0.1, "+x0:%.5f" % (0.030 * s), "TKSUP", i, mirror=True)
    tec = [127.75, 141.75, 155.75, 169.75, 183.75, 202.25, 222.75, 245.25, 269.75]
    for i, z in enumerate(tec, 1):
        tube("TEC%d" % i, "G4_Si", 22.0, 113.5, z, "+um:400", "TEC", i, sens=True, mirror=True)
        tube("TEC%dS" % i, "TkSupport", 22.0, 113.5, z + 0.1, "+x0:%.5f" % (0.040 * s), "TKSUP", i, mirror=True)
    # services and end flanges: the material peak near |eta| = 1.4 [M]
    tube("TIBSERV", "TkServ", 20.0, 56.0, 70.5, "+x0:%.5f" % (0.60 * s), "TKSUP", mirror=True)
    tube("TOBSERV", "TkSupport", 57.0, 113.5, 111.0, "+x0:%.5f" % (0.30 * s), "TKSUP", mirror=True)

# ---------------- ECAL [D]: PbWO4, EB r 129-152 cm (25.8 X0), EE 22 cm (24.7 X0), preshower
w("# ECAL [D]")
c = a.calo_scale
tube("EB", "G4_PbWO4", 129.0, "+cm:%.3f" % (23.0 * c), -290, 290, "EB", sens=True, calo=True)
if not a.no_endcaps:
    tube("ESPB1", "G4_Pb", 45.0, 123.0, 298.0, "+x0:2", "ES", 1, calo=True, mirror=True)
    tube("ESSI1", "G4_Si", 45.0, 123.0, 299.3, "+um:320", "ES", 1, sens=True, calo=True, mirror=True)
    tube("ESPB2", "G4_Pb", 45.0, 123.0, 299.5, "+x0:1", "ES", 2, calo=True, mirror=True)
    tube("ESSI2", "G4_Si", 45.0, 123.0, 300.2, "+um:320", "ES", 2, sens=True, calo=True, mirror=True)
    tube("EE", "G4_PbWO4", 32.0, 171.0, 317.0, "+cm:%.3f" % (22.0 * c), "EE", sens=True, calo=True, mirror=True)

# ---------------- HCAL [D]: brass/scintillator sampling, HB r 177-295 cm (5.8 lambda_I at eta = 0)
w("# HCAL [D]: 3.7 mm scintillator layers; HB brass 15 x 5.5 cm between steel plates of 4 and 7.5 cm")
r = 177.0
tube("HBF", "G4_STAINLESS-STEEL", r, "+cm:%.3f" % (4.0 * c), -390, 390, "HB", calo=True)
r += 4.0 * c + 0.1
for i in range(16):
    tube("HBS%d" % i, "G4_PLASTIC_SC_VINYLTOLUENE", r, "+cm:0.37", -390, 390, "HB", i + 1, sens=True, calo=True)
    r += 0.37 + 0.53
    if i < 15:
        tube("HBA%d" % i, "G4_BRASS", r, "+cm:%.3f" % (5.5 * c), -390, 390, "HB", calo=True)
        r += 5.5 * c + 0.05
tube("HBB", "G4_STAINLESS-STEEL", r, "+cm:%.3f" % (7.5 * c), -390, 390, "HB", calo=True)
assert r + 7.5 * c < 295.0, r
if not a.no_endcaps:
    z = 391.0
    for i in range(17):
        tube("HES%d" % i, "G4_PLASTIC_SC_VINYLTOLUENE", 40.0, 295.0, z, "+cm:0.37", "HE", i + 1, sens=True,
             calo=True, mirror=True)
        z += 0.37 + 0.23
        tube("HEA%d" % i, "G4_BRASS", 40.0, 295.0, z, "+cm:%.3f" % (7.9 * c), "HE", calo=True, mirror=True)
        z += 7.9 * c + 0.1
    assert z < 575.0, z

# ---------------- solenoid [D]: cryostat walls + 31 cm cold mass, 12.5 m long (about 1.1 lambda_I)
w("# solenoid [D]")
tube("CRYO_IN", "G4_STAINLESS-STEEL", 297.0, "+cm:2", -625, 625, "COIL")
tube("COIL", "CoilMix", 318.0, "+cm:31", -625, 625, "COIL")
tube("CRYO_OUT", "G4_STAINLESS-STEEL", 375.0, "+cm:3", -625, 625, "COIL")

# ---------------- HO and the barrel return yoke with 4 DT stations (approximate radii) [D]
w("# HO, barrel yoke (iron 30 + 63 + 63 cm) and muon stations MB1-4, two gas layers each; approximate radii")
tube("HO0", "G4_PLASTIC_SC_VINYLTOLUENE", 382.0, "+cm:1", -126, 126, "HO", 1, sens=True, calo=True)
tube("TC", "G4_Fe", 384.0, "+cm:19.5", -126, 126, "HO", calo=True)
tube("HO1", "G4_PLASTIC_SC_VINYLTOLUENE", 407.0, "+cm:1", -635, 635, "HO", 2, sens=True, calo=True)
yoke = [(460.0, 30.0), (540.0, 63.0), (650.0, 63.0)]
mb = [420.0, 500.0, 610.0, 720.0]
for i, rr in enumerate(mb, 1):
    for k in range(2):
        tube("MB%d%d" % (i, k), "MuonGas", rr + 20.0 * k, "+cm:2", -635, 635, "MB", i, sens=True)
for i, (rr, t) in enumerate(yoke, 1):
    tube("YB%d" % i, "G4_Fe", rr, "+cm:%.1f" % t, -635, 635, "YOKE", i)
if not a.no_endcaps:
    w("# endcap muon stations ME1-4 and yoke disks YE1-3 (60, 60, 25 cm, r > 1 m); approximate positions")
    tube("ME11", "MuonGas", 100.0, 275.0, 585.0, "+cm:2", "ME", 1, sens=True, mirror=True)
    tube("ME12", "MuonGas", 275.0, 700.0, 690.0, "+cm:2", "ME", 1, sens=True, mirror=True)
    tube("YE1", "G4_Fe", 100.0, 700.0, 730.0, "+cm:60", "YOKE", 4, mirror=True)
    tube("ME2", "MuonGas", 140.0, 700.0, 820.0, "+cm:2", "ME", 2, sens=True, mirror=True)
    tube("YE2", "G4_Fe", 100.0, 700.0, 840.0, "+cm:60", "YOKE", 5, mirror=True)
    tube("ME3", "MuonGas", 160.0, 700.0, 920.0, "+cm:2", "ME", 3, sens=True, mirror=True)
    tube("YE3", "G4_Fe", 100.0, 700.0, 945.0, "+cm:25", "YOKE", 6, mirror=True)
    tube("ME4", "MuonGas", 160.0, 700.0, 1010.0, "+cm:2", "ME", 4, sens=True, mirror=True)

# ---------------- field: uniform regions in (r, |z|), first match wins
w("# field: Bz[T] rmin rmax |z|min |z|max")
if a.field != "none":
    w("field 3.8 0 315 0 600")
if a.field == "yoke":
    for rr, t in yoke:
        w("field -1.8 %.1f %.1f 0 635" % (rr, rr + t))

# ---------------- boundaries: first exit time of each quirk is written out
w("# boundaries: name r |z|")
for b in [("PIX", 17, 52), ("TRK", 115, 280), ("ECALIN", 129, 315), ("HCALIN", 177, 390), ("COIL", 297, 600),
          ("MUON", 420, 690), ("YOKEOUT", 715, 1000)]:
    w("boundary %s %d %d" % b)

open(a.out, "w").write("\n".join(L) + "\n")
print("wrote", a.out, tag)
