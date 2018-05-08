#! /bin/bash

##### Initialize Benchmarks #####
source xeon_scripts/common_variables.sh
make distclean

##### Launch Tests #####
echo "Tar and send to KNL"
./xeon_scripts/tarAndSendToKNL.sh

echo "Run benchmarking on KNL concurrently with SKL-SP benchmarks" 
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build-knl.sh >& benchmark_knl_dump.txt &

echo "Tar and send to SNB"
./xeon_scripts/tarAndSendToSNB.sh

echo "Run benchmarking on SNB concurrently with SKL-SP benchmarks" 
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build-snb.sh >& benchmark_snb_dump.txt &

echo "Run benchmarking on SKL-SP"
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build.sh SKL-SP

##### Validation tests #####
echo "Running ROOT based validation"
./val_scripts/validation-cmssw-benchmarks.sh

echo "Waiting for KNL and SNB"
wait

##### Benchmark Plots #####
echo "Producing benchmarking plots"
./plotting/benchmarkPlots.sh

##### Plots from Text Files #####
echo "Producing plots from text files"
./plotting/textDumpPlots.sh

##### Final cleanup #####
make distclean
