#!/bin/bash
# standalone Combine on LCG_106 (no CMSSW)
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh
cd $(dirname $0)/HiggsAnalysis/CombinedLimit
rm -rf build; mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_VDT=OFF .. && make -j3 && make install
