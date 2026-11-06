#!/bin/bash
# run_point.sh <model taup|quirk|quirkb|quirks|quirksb> <mass GeV> <lambda eV> <ngen> <seed> <outdir> : gen, filter, quirk-only sim
M=$1; MASS=$2; L=$3; N=$4; SEED=$5; OUT=$6
mkdir -p $OUT
R=$(cd $(dirname $0)/..; pwd)
B=${QP_BIN:-$R/build}  # quirksim build
C=${QP_CARDS:-$R/cards}
G=${QP_GEN:-$R/gen}
export LD_LIBRARY_PATH=$B/quirkcore:$LD_LIBRARY_PATH
T=${M}_m${MASS}_L${L}_s${SEED}
W=$(mktemp -d ${TMPDIR:-/tmp}/qp_XXXX)
cd $W
CMND=$G/pythia_quirk_dy_13TeV.cmnd  # SU(2) doublet
case $M in taup|quirks*) CMND=$G/pythia_taup1e_13TeV.cmnd;; esac  # SU(2) singlet (EXO-18-002 tau' card)
$B/quirksim/quirkgen $MASS $N gen.hepmc3 $CMND $SEED 35 > $OUT/$T.gen.log 2>&1 || exit 1
F=""; [[ $M == quirk*b ]] && F="1.3 bound"  # bound-pair GEN filter
python3 $G/prep_gen.py gen.hepmc3 genf.hepmc3 $OUT/$T.npz $F >> $OUT/$T.gen.log 2>&1 || exit 2
NK=$(grep -c "^E " genf.hepmc3)
cat > sim.mac <<EOM
/quirk/mass $MASS GeV
/quirk/lambda $L eV
/qsim/card $C/cms_like.card
/qsim/gen/mode hepmc
/qsim/gen/hepmcFile genf.hepmc3
/qsim/gen/quirksOnly true
/qsim/out/towers true
/qsim/out/printEvery 50
/run/initialize
/run/beamOn $NK
EOM
$B/quirksim/quirksim -m sim.mac -t ${QP_THREADS:-1} -s $SEED -o $OUT/$T.root > $OUT/$T.sim.log 2>&1 || exit 3
cd /; rm -rf $W
echo done > $OUT/$T.done
