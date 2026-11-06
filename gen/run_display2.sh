#!/bin/bash
# display gun pair (m 200 GeV) through the CMS-like detector, all steps
S=$(cd $(dirname $0)/..; pwd); O=${QP_OUT:-$S/out}/display; mkdir -p $O
export LD_LIBRARY_PATH=$S/build/quirkcore:$LD_LIBRARY_PATH
for L in 30 100 300; do
cat > $O/ed_L$L.mac <<EOM
/quirk/mass 200 GeV
/quirk/lambda $L eV
/qsim/card $S/cards/cms_like.card
/qsim/gen/mode gun
/qsim/gun/p1 300 50 100 GeV
/qsim/gun/p2 250 -40 150 GeV
/qsim/gun/vertex 0 0 0 mm
/qsim/out/stepPrescale 1
/run/initialize
/run/beamOn 1
EOM
$S/build/quirksim/quirksim -m $O/ed_L$L.mac -t 1 -s 3 -o $O/ed_L$L.root > $O/ed_L$L.log 2>&1 &
done; wait
