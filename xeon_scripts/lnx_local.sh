#! /bin/bash

./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build.sh LNX-G forPR
./val_scripts/validation-cmssw-benchmarks.sh forPR
./plotting/benchmarkPlots_lnx.sh forPR
make distclean
