#!/bin/bash
# run_limits.sh <model> <m1,m2,..> <L1,L2,..> [tag] : cards at sigma_LO and asymptotic limits (r = sigma95 / sigma_LO)
source $(cd $(dirname $0)/..; pwd)/ext/env_combine.sh
MODEL=$1; MASSES=${2//,/ }; LAMS=${3//,/ }; TAG=${4:-nom}  # comma-separated lists
cd $(dirname $0); OUTDIR=${QP_OUT:-$(cd ..; pwd)/out}; CARDDIR=$OUTDIR/cards/$TAG; mkdir -p $CARDDIR
for mass in $MASSES; do for lam in $LAMS; do
  name=${MODEL}_m${mass}_L${lam}; emu=$OUTDIR/emu/$name.npz; [ $TAG != nom ] && emu=$OUTDIR/emu/$TAG/$name.npz
  [ -f $emu ] || continue
  [ -f $CARDDIR/lim_$name.log ] && [ $CARDDIR/lim_$name.log -nt $emu ] && continue  # up to date
  simdir=$OUTDIR/$MODEL; [ $MODEL != taup ] && simdir=$OUTDIR/$MODEL/flat
  sigma=$(python3 xsec.py $simdir "${MODEL}_m${mass}_L${lam}_s*")
  python3 make_card.py $emu $name $sigma $CARDDIR > $CARDDIR/card_$name.mklog 2>&1 || continue
  (cd $CARDDIR; combine -M AsymptoticLimits card_$name.txt -n $name > lim_$name.log 2>&1
   exp_hi=$(awk '/Expected 97.5%/{print $NF}' lim_$name.log); obs=$(awk '/Observed Limit/{print $NF}' lim_$name.log)
   if [ -n "$exp_hi" ] && [ -n "$obs" ] && awk "BEGIN{exit !($obs > 2 * $exp_hi)}"; then  # failed observed fit: retry with a bounded r range
     combine -M AsymptoticLimits card_$name.txt -n $name --rMax $(awk "BEGIN{print 10 * $exp_hi}") > lim_$name.log 2>&1
   fi
   echo "$mass $lam $sigma" > xs_$name.txt)
done; done
