#!/bin/bash
# MIP reference: charge-1 10 GeV "quirks" (Lambda 1 eV, free) from a gun, p = 700 GeV (beta gamma 70), several eta
S=$(cd $(dirname $0)/..; pwd); O=/ceph/cms/store/user/tvami/HSCP-Quirks/QuirkPheno/mip
export LD_LIBRARY_PATH=$S/build/quirkcore:$LD_LIBRARY_PATH
for eta in 0.05 0.3 0.6 0.9; do
read a b c d e f <<< $(python3 -c "
import math; p=700.; e=$eta; t=2*math.atan(math.exp(-e)); s=math.sin(t); z=math.cos(t)
print(p*s*math.cos(0.4), p*s*math.sin(0.4), p*z, p*s*math.cos(2.9), p*s*math.sin(2.9), -p*z)")
cat > $O/mip_eta$eta.mac <<EOM
/quirk/mass 10 GeV
/quirk/lambda 1 eV
/qsim/card $S/cards/cms_like.card
/qsim/gen/mode gun
/qsim/gun/p1 $a $b $c GeV
/qsim/gun/p2 $d $e $f GeV
/qsim/gun/vertex 0 0 0 mm
/qsim/out/printEvery 500
/run/initialize
/run/beamOn ${NMIP:-2000}
EOM
$S/build/quirksim/quirksim -m $O/mip_eta$eta.mac -t 2 -s 7 -o $O/mip_eta$eta.root > $O/mip_eta$eta.log 2>&1 &
done; wait
