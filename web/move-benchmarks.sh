#! /bin/bash

# command line input
dir=${1:-"benchmarks"} # Main output dir name
suite=${2:-"forPR"} # which set of benchmarks to run: full, forPR, forConf
afs_or_eos=${3:-"afs"} # which user space to use: afs or eos

# source global variables
source xeon_scripts/common-variables.sh ${suite}
source xeon_scripts/init-env.sh

# First collect all plots and text files into common dir
echo "Moving plots and text files locally to ${dir}"
./web/collectBenchmarks.sh ${dir} ${suite}

# Now copy to lxplus
echo "Moving plots and text files remotely to lxplus"
./web/tarAndSendToLXPLUS.sh ${dir} ${suite} ${afs_or_eos}

# Final cleanup of directory
echo "Removing local files"
./xeon_scripts/trashSKL-SP.sh 
rm -rf ${dir}
