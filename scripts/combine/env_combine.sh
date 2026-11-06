# source me: standalone Combine (LCG_106)
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh
I=$(cd $(dirname ${BASH_SOURCE[0]}); pwd)/HiggsAnalysis/CombinedLimit/install
export PATH=$I/bin:$PATH LD_LIBRARY_PATH=$I/lib:$LD_LIBRARY_PATH PYTHONPATH=$I/lib/python3.11/site-packages:$PYTHONPATH
