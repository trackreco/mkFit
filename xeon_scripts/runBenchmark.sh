#! /bin/bash

##### Initialize Benchmarks #####
source xeon_scripts/common_variables.sh
make distclean

##### Launch Tests #####
echo "Tar and send to KNL"
./xeon_scripts/tarAndSendToKNL.sh

echo "Run benchmarking on KNL concurrently with SNB benchmarks" 
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build-knl.sh >& benchmark_knl_dump.txt &

echo "Run benchmarking on SNB"
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build.sh SNB

echo "Waiting for KNL"
wait

##### Benchmark Plots #####
./plotting/benchmarkPlots.sh

##### Plots from Text Files #####
./plotting/textDumpPlots.sh

##### Validation tests #####
echo "Running ROOT based validation"
./val_scripts/validation-snb-cmssw-ttbar-fulldet-build-allval.sh

##### Final cleanup #####
make distclean
