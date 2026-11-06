# LCG_110 (Geant4 11.4.1, HepMC3, Pythia8, ROOT) from /cvmfs/sft.cern.ch; el9 only (scripts/el9 elsewhere)
source /cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-el9-gcc13-opt/setup.sh
export PYTHONPATH=${QP_PYLIB:-$(cd $(dirname ${BASH_SOURCE[0]})/..; pwd)/ext/pylib}:$PYTHONPATH  # cmsstyle (pip --target)
export HDF5_USE_FILE_LOCKING=FALSE
