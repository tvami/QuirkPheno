#!/bin/bash
# tau' 1e limit closure: emulated signal at the paper's sigma_ref on the published background model
source /home/users/tvami/HSCP-Quirks/QuirkPheno/ext/env_combine.sh
cd $(dirname $0); B=/ceph/cms/store/user/tvami/HSCP-Quirks/QuirkPheno; O=$B/cards/closure; mkdir -p $O
for m in 800 1000 1400 1800 2200 2600; do
  s=$(python3 -c "import json; print(json.load(open('card_yields.json'))['$m']['sigma_ref'])")
  python3 make_card.py $B/emu/taup_m${m}_L1.npz taup$m $s $O
  (cd $O; combine -M AsymptoticLimits card_taup$m.txt -n taup$m > lim_taup$m.log 2>&1)
done
