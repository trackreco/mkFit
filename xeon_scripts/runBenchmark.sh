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
echo "Producing benchmarking plots"
./plotting/benchmarkPlots.sh

##### Plots from Text Files #####
echo "Producing plots from text files"
./plotting/textDumpPlots.sh

##### Validation tests #####
echo "Running ROOT based validation"
./val_scripts/validation-snb-cmssw-benchmarks.sh

##### Final cleanup #####
make distclean
