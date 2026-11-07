#!/bin/bash
# emu_all.sh <model> <mass list> [lambda list] : emulate every point (patterns over seeds)
MODEL=$1; MASSES=$2; LAMS=${3:-1}
OUTDIR=${QP_OUT:-$(cd $(dirname $0)/..; pwd)/out}; cd $(dirname $0)
for mass in $MASSES; do for lam in $LAMS; do
  simdir=$OUTDIR/$MODEL; pattern="${MODEL}_m${mass}_L${lam}_s*"
  [ $MODEL != taup ] && simdir=$simdir/flat  # linked by finalize.sh
  emu=$OUTDIR/emu/${MODEL}_m${mass}_L${lam}.npz
  if [ -f $emu ] && [ -z "$(find -L $simdir -name "${pattern}.root" -newer $emu)" ]; then continue; fi  # up to date
  python3 emulate.py $simdir "$pattern" $emu
done; done
