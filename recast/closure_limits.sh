#!/bin/bash
# tau' 1e limit closure: emulated signal at the paper's sigma_ref on the published background model
source $(cd $(dirname $0)/..; pwd)/ext/env_combine.sh
cd $(dirname $0); OUTDIR=${QP_OUT:-$(cd ..; pwd)/out}; CARDDIR=$OUTDIR/cards/closure; mkdir -p $CARDDIR
for mass in 800 1000 1400 1800 2200; do
  sigma=$(python3 -c "import json; print(json.load(open('card_yields.json'))['$mass']['sigma_ref'])")
  python3 make_card.py $OUTDIR/emu/taup_m${mass}_L1.npz taup$mass $sigma $CARDDIR
  (cd $CARDDIR; combine -M AsymptoticLimits card_taup$mass.txt -n taup$mass > lim_taup$mass.log 2>&1)
done
