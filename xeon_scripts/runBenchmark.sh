#! /bin/bash

##### Command Line Input #####
suite=${1:-"forPR"} # which set of benchmarks to run: full, forPR, forConf
useLNX=${2:-0}
lnxuser=${3:-${USER}}

##### Initialize Benchmarks #####
source xeon_scripts/common-variables.sh ${suite} ${useLNX} ${lnxuser}
source xeon_scripts/init-env.sh
make distclean

##### Launch Tests #####
if [[ ${useLNX} -eq 1 ]] || [[ ${useLNX} -eq 2 ]]
then
echo "Tar and send to LNX7188"
./xeon_scripts/tarAndSendToRemote.sh LNX-G ${suite} ${useLNX} ${lnxuser}

echo "Run benchmarking on LNX7188 concurrently with SKL-SP benchmarks" 
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build-remote.sh LNX-G ${suite} ${useLNX} ${lnxuser} >& benchmark_lnx-g_dump.txt &

echo "Tar and send to LNX4108"
./xeon_scripts/tarAndSendToRemote.sh LNX-S ${suite} ${useLNX} ${lnxuser}

echo "Run benchmarking on LNX4108 concurrently with SKL-SP benchmarks" 
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build-remote.sh LNX-S ${suite} ${useLNX} ${lnxuser} >& benchmark_lnx-s_dump.txt &
fi

if [[ ${useLNX} -eq 0 ]] || [[ ${useLNX} -eq 2 ]]
then

echo "Tar and send to KNL"
./xeon_scripts/tarAndSendToRemote.sh KNL ${suite} ${useLNX} ${lnxuser}

echo "Run benchmarking on KNL concurrently with SKL-SP benchmarks" 
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build-remote.sh KNL ${suite} ${useLNX} ${lnxuser} >& benchmark_knl_dump.txt &

echo "Tar and send to SNB"
./xeon_scripts/tarAndSendToRemote.sh SNB ${suite} ${useLNX} ${lnxuser}

echo "Run benchmarking on SNB concurrently with SKL-SP benchmarks" 
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build-remote.sh SNB ${suite} ${useLNX} ${lnxuser} >& benchmark_snb_dump.txt &

echo "Run benchmarking on SKL-SP"
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build.sh SKL-SP ${suite} ${useLNX} ${lnxuser} 
fi

####### Validation tests #####
echo "Running ROOT based validation"
./val_scripts/validation-cmssw-benchmarks.sh ${suite}

echo "Waiting for LNX-G, LNX-S, KNL and SNB"
wait

##### Benchmark Plots #####
echo "Producing benchmarking plots"
./plotting/benchmarkPlots.sh ${suite} ${useLNX} ${lnxuser} 

##### Plots from Text Files #####
echo "Producing plots from text files"
./plotting/textDumpPlots.sh ${suite} ${useLNX} ${lnxuser} 

##### Final cleanup #####
make distclean

##### Final message #####
echo "Finished benchmarking and validation suite!"
