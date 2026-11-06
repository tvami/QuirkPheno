# Detector cards

A card is a list of cylinders and disks placed directly in a cylindrical world (no hierarchy,
overlaps are checked at construction). Lengths in cm, fields in tesla, `#` starts a comment.

| Keyword | Syntax | Meaning |
|---|---|---|
| `world` | `world <r> <half z> <material>` | must come first |
| `checkOverlaps` | `checkOverlaps 0/1` | Geant4 overlap check of every placement (default 1) |
| `mix` | `mix <name> <density g/cm3> <material> <mass fraction> ...` | mixture of NIST or earlier materials |
| `scaled` | `scaled <name> <material> <density factor>` | material with a scaled density |
| `tube` | `tube <name> <material> <rmin> <rmax> <zmin> <zmax> [det=<group>] [layer=<n>] [sens] [calo] [mirror]` | volume; `rmax`/`zmax` may be a thickness from the inner edge: `+x0:<X/X0>`, `+lam:<lambda_I>`, `+cm:<cm>`, `+um:<um>`; `mirror` adds the copy at negative z |
| `field` | `field <Bz> <rmin> <rmax> <abs z min> <abs z max>` | uniform Bz in a region; the first matching region wins, zero elsewhere |
| `boundary` | `boundary <name> <r> <half z>` | the first time each quirk leaves this cylinder is written out |

`sens` volumes record quirk hits; `calo` volumes add their energy to the 0.087 x 0.087 towers;
`det=` groups volumes for the per-group energies (default: the volume name).

## Cards

`cms_like.card` is the generic CMS-like detector, written by `make_cms_like.py`.

## The CMS-like card: sources and approximations

All numbers are from public documents, rounded:

- [D] CMS Collaboration, "The CMS experiment at the CERN LHC", JINST 3 (2008) S08004;
- [P] CMS Collaboration, "CMS Technical Design Report for the Pixel Detector Upgrade", CERN-LHCC-2012-016, CMS-TDR-011.

| System | Model | Source |
|---|---|---|
| Beam pipe | Be, r = 2.17-2.25 cm, \|z\| < 3 m, then 1.5 mm steel to \|z\| = 11 m (approximate, the real pipe is conical); vacuum inside | [D, P] |
| Pixel | 4 barrel layers at r = 2.9, 6.8, 10.9, 16.0 cm (\|z\| < 26.7 cm), 3 disks per side at \|z\| = 29.1, 39.6, 51.6 cm; 285 um Si + support | [P] |
| Strips | TIB r = 25.5, 33.9, 41.85, 49.8 cm (\|z\| < 70 cm); TOB r = 60.8, 69.2, 78.0, 86.8, 96.5, 108.0 cm (\|z\| < 109 cm); layers 1-2 double sided; TID 3 and TEC 9 disks per side | [D] |
| Tracker material | support layers of a C/Al/Cu/Kapton mixture, X/X0 set per layer, plus service disks; total 0.374 X0 at eta = 0 and a peak of 1.78 X0 at \|eta\| = 1.1 (public: about 0.4 and 1.8 near 1.4, [D] Fig. 3.4) | tuned to [D] |
| ECAL | PbWO4, EB r = 129-152 cm (25.8 X0, \|z\| < 2.9 m), EE 22 cm (24.7 X0) at \|z\| = 3.17 m, preshower Pb 2 + 1 X0 with Si | [D] |
| HCAL | HB r = 177-295 cm: steel 4 cm, 15 brass plates of 5.5 cm, 16 scintillator layers of 3.7 mm, steel 7.5 cm (5.6 lambda_I with the scintillator); HE 17 brass 7.9 cm plates from \|z\| = 3.91 m | [D] |
| Solenoid | cryostat walls of steel (2 and 3 cm) and a 31 cm Al-dominated cold mass, r = 2.97-3.78 m, \|z\| < 6.25 m (about 1.1 lambda_I) | [D] |
| HO | scintillator at r = 3.82 m (ring 0) and 4.07 m, tail catcher 19.5 cm Fe in ring 0 | [D] |
| Muon barrel | 3 iron layers (30, 63, 63 cm) at r = 4.6, 5.4, 6.5 m, 4 stations of 2 Ar/CO2 gas layers at r = 4.2, 5.0, 6.1, 7.2 m, \|z\| < 6.35 m | approximate |
| Muon endcap | yoke disks of 60, 60, 25 cm for r > 1 m, 4 gas stations between \|z\| = 5.85 and 10.1 m | approximate |
| Field | 3.8 T for r < 3.15 m and \|z\| < 6 m; -1.8 T in the barrel iron; zero elsewhere (endcap return field not modeled) | simplified |

The radii and z positions of the muon system and the endcap yoke are rounded estimates, not
taken from a table: check them against [D] (Figs. 1.2, 7.x) before any of them enters a paper.
The calorimeter absorber is modeled as homogeneous plates; the energy written in `eDet` is all
deposited energy, `eDetSens` only the energy in the active layers (no Birks law, no sampling
calibration).
