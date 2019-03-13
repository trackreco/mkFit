#! /bin/bash

##### Command Line Input #####
suite=${1:-"forPR"} # which set of benchmarks to run: full, forPR, forConf

##### Initialize Benchmarks #####
source xeon_scripts/common-variables.sh ${suite}
source xeon_scripts/init-env.sh
make distclean

##### Launch Tests #####
echo "Tar and send to KNL"
./xeon_scripts/tarAndSendToRemote.sh KNL ${suite}

echo "Run benchmarking on KNL concurrently with SKL-SP benchmarks" 
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build-remote.sh KNL ${suite} >& benchmark_knl_dump.txt &

echo "Tar and send to SNB"
./xeon_scripts/tarAndSendToRemote.sh SNB ${suite}

echo "Run benchmarking on SNB concurrently with SKL-SP benchmarks" 
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build-remote.sh SNB ${suite} >& benchmark_snb_dump.txt &

echo "Run benchmarking on SKL-SP"
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build.sh SKL-SP ${suite}

##### Validation tests #####
echo "Running ROOT based validation"
./val_scripts/validation-cmssw-benchmarks.sh ${suite}

echo "Waiting for KNL and SNB"
wait

##### Benchmark Plots #####
echo "Producing benchmarking plots"
./plotting/benchmarkPlots.sh ${suite}

##### Plots from Text Files #####
echo "Producing plots from text files"
./plotting/textDumpPlots.sh ${suite}

##### Final cleanup #####
make distclean

##### Final message #####
echo "Finished benchmarking and validation suite!"
