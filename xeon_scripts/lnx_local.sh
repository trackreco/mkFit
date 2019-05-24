#! /bin/bash

suite=${1:-"forPR"}
lnxuser=${2:-${USER}}
useLNX=${3:-1} #0=phi only, 1=lnx only, 2= phi and lnx

##### Initialize Benchmarks #####
source xeon_scripts/common-variables.sh ${suite} ${lnxuser}
source xeon_scripts/init-env.sh
make distclean

##### Start Tests######
echo "Tar and send to LNX4108"
./xeon_scripts/tarAndSendToRemote.sh LNX-S ${suite} ${lnxuser}

echo "Run benchmarking on LNX4108 concurrently with SKL-SP benchmarks"
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build-remote.sh LNX-S ${suite} ${lnxuser} >& benchmark_lnx-s_dump.txt &

echo "Run Benchmarking on LNX7188"
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build.sh LNX-G forPR

####### Validation test on LNX7188#######
#./val_scripts/validation-cmssw-benchmarks.sh forPR
#echo "waiting for LNX4108"
wait

########Plotting
echo "producing benchmarking pltos"
./plotting/benchmarkPlots_lnx.sh forPR ${suite} ${lnxuser} ${useLNX}

###### Plots from Text Files #####
echo "Producing plots from text files"
./plotting/textDumpPlots_lnx.sh ${suite} ${lnxuser} ${useLNX}

####final cleaning####
make distclean
echo "finished benchmarking and validation" 
