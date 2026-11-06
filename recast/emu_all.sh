#!/bin/bash
# emu_all.sh <model> <mass list> [lambda list] : emulate every point (patterns over seeds)
M=$1; MASSES=$2; LAMS=${3:-1}
B=/ceph/cms/store/user/tvami/HSCP-Quirks/QuirkPheno; cd $(dirname $0)
for m in $MASSES; do for L in $LAMS; do
  D=$B/$M; P="${M}_m${m}_L${L}_s*"
  [ $M != taup ] && D=$D/flat  # linked by finalize.sh
  O=$B/emu/${M}_m${m}_L${L}.npz
  if [ -f $O ] && [ -z "$(find -L $D -name "${P}.root" -newer $O)" ]; then continue; fi  # up to date
  python3 emulate.py $D "$P" $O
done; done
