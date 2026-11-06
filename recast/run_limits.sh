#!/bin/bash
# run_limits.sh <model> <m1,m2,..> <L1,L2,..> [tag] : cards at sigma_LO and asymptotic limits (r = sigma95 / sigma_LO)
source /home/users/tvami/HSCP-Quirks/QuirkPheno/ext/env_combine.sh
M=$1; MASSES=${2//,/ }; LAMS=${3//,/ }; TAG=${4:-nom}  # comma-separated lists
cd $(dirname $0); B=/ceph/cms/store/user/tvami/HSCP-Quirks/QuirkPheno; O=$B/cards/$TAG; mkdir -p $O
for m in $MASSES; do for L in $LAMS; do
  n=${M}_m${m}_L${L}; E=$B/emu/$n.npz; [ $TAG != nom ] && E=$B/emu/$TAG/$n.npz
  [ -f $E ] || continue
  [ -f $O/lim_$n.log ] && [ $O/lim_$n.log -nt $E ] && continue  # up to date
  D=$B/$M; [ $M != taup ] && D=$B/$M/flat
  s=$(python3 xsec.py $D "${M}_m${m}_L${L}_s*")
  python3 make_card.py $E $n $s $O > $O/card_$n.mklog 2>&1 || continue
  (cd $O; combine -M AsymptoticLimits card_$n.txt -n $n > lim_$n.log 2>&1
   e=$(awk '/Expected 97.5%/{print $NF}' lim_$n.log); o=$(awk '/Observed Limit/{print $NF}' lim_$n.log)
   if [ -n "$e" ] && [ -n "$o" ] && awk "BEGIN{exit !($o > 2 * $e)}"; then  # failed observed fit: retry with a bounded r range
     combine -M AsymptoticLimits card_$n.txt -n $n --rMax $(awk "BEGIN{print 10 * $e}") > lim_$n.log 2>&1
   fi
   echo "$m $L $s" > xs_$n.txt)
done; done
