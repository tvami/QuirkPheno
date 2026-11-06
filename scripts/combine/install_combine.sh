#!/bin/bash
# install step with a local python site dir
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh
cd $(dirname $0)/HiggsAnalysis/CombinedLimit/build
I=$(cd .. && pwd)/install
cmake -DPython_SITELIB=$I/lib/python3.11/site-packages .. > /dev/null && make install 2>&1 | tail -2
