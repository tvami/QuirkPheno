Standalone Combine v10.2.1 on LCG_106 (el9). Clone into `ext/HiggsAnalysis/CombinedLimit`, copy these scripts to `ext/`,
run `build_combine.sh` (`scripts/el9 -b` on non-el9 hosts), then `install_combine.sh`. ROOT 6.32 fix in `python/ShapeTools.py`:
`histpdf.get().second()` -> `histpdf.get()[1]`.
