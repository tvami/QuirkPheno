# QuirkPheno

Code for "Constraints on quirks from a Geant4-based reinterpretation of the CMS search for heavy stable charged
particles": a reinterpretation of the CMS search EXO-18-002 (arXiv:2410.09164, HEPData ins2840007) for Drell-Yan
produced quirks, using only public information.

## Requirements

- CVMFS with `sft.cern.ch`: LCG_110 (Geant4 11.4.1, HepMC3, Pythia8, ROOT) and LCG_106 (Combine).
- An el9 host, or apptainer: `scripts/el9 <command>` runs in the public image
  `/cvmfs/singularity.opensciencegrid.org/htc/rocky:9` with LCG_110 set up.
- `cmsstyle` for the plots: `pip install --target ext/pylib cmsstyle`.

## Build

```
scripts/el9 "cmake -S . -B build && cmake --build build -j4"   # quirkcore, quirksim, quirkgen
```

Combine: see `scripts/combine/README.md`.

## Chain

Outputs go to `$QP_OUT` (default `out/`).

```
hepdata/fetch.sh                                              # published likelihoods and efficiency maps
gen/run_mip.sh; recast/build_mip.py                           # MIP dE/dx templates
gen/run_point.sh taup <m> 1 <n> <seed> $QP_OUT/taup           # tau' closure sample
gen/run_point.sh quirk <m> <L> <n> <seed> $QP_OUT/quirk/flat  # quirkb: bound-pair filter; quirks, quirksb: singlet
recast/emu_all.sh <model> "<masses>" "<lambdas>"              # emulated reconstruction and selection
recast/card_yields.py; recast/closure_limits.sh               # tau' closure
recast/run_limits.sh <model> <m1,m2,..> <L1,L2,..>            # quirk limits
recast/collect.py nom doublet                                 # or: nom singlet
plots/limits.py; plots/faser_cmp.py; plots/closure.py         # figures in plots/figs
```

| Directory | Content |
|---|---|
| `quirkcore/`, `quirksim/` | quirk transport library and the Geant4 application |
| `cards/` | CMS-like detector card and its generator |
| `gen/` | Pythia8 cards, generation, generator-level filter, quirk-only simulation |
| `recast/` | emulated reconstruction and selection, datacards, limits |
| `plots/` | figure scripts |
| `hepdata/` | download of the HEPData record, FASER contours |
| `scripts/` | el9 wrapper, LCG environment, Combine build |

## License

Apache-2.0 (`LICENSE`); see `NOTICE` for the code `quirkcore` is derived from.
