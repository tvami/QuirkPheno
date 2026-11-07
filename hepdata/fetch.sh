#!/bin/bash
# fetch the original HEPData submission of EXO-18-002 (yaml + resources)
cd $(dirname $0)
curl -sL -o original.zip "https://www.hepdata.net/download/submission/ins2840007/1/original"
mkdir -p original && (cd original && unzip -qo ../original.zip)
