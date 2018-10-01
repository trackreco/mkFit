#! /bin/bash

###########
## Input ##
###########

ben_arch=SKL-SP

###################
## Configuration ##
###################

## Source environment and common variables
source xeon_scripts/common-variables.sh ${suite}
source xeon_scripts/init-env.sh

## Platform specific settings
if [[ "${ben_arch}" == "SNB" ]]
then
    mOpt="-j 12"
    dir=/data2/nfsmic/slava77/samples
    maxth=24
    maxvu=8
    declare -a nths=("1" "2" "4" "6" "8" "12" "16" "20" "24")
    declare -a nvus=("1" "2" "4" "8")
    declare -a nevs=("1" "2" "4" "8" "12")
elif [[ "${ben_arch}" == "KNL" ]]
then
    mOpt="-j 64 AVX_512:=1"
    dir=/data1/work/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000
    maxth=256
    maxvu=16
    declare -a nths=("1" "2" "4" "8" "16" "32" "64" "96" "128" "160" "192" "224" "256")
    declare -a nvus=("1" "2" "4" "8" "16")
    declare -a nevs=("1" "2" "4" "8" "16" "32" "64" "128")
elif [[ "${ben_arch}" == "SKL-SP" ]]
then
    mOpt="-j 32 AVX_512:=1"
    dir=/data2/slava77/samples
    maxth=64
    maxvu=16
    declare -a nths=("1" "2" "4" "8" "16" "32" "48" "64")
    declare -a nvus=("1" "2" "4" "8" "16")
    declare -a nevs=("1" "2" "4" "8" "16" "32" "64")
else 
    echo ${ben_arch} "is not a valid architecture! Exiting..."
    exit
fi

## Common file setup
subdir=2017/pass-4874f28/initialStep/10muEtaLT06Pt1to10/
file=memoryFile.fv3.clean.writeAll.recT.072617.bin
nevents=500

## Common executable setup
minth=1
minvu=1
seeds="--cmssw-n2seeds"
exe="./mkFit/mkFit ${seeds} --input-file ${dir}/${subdir}/${file}"

## Common output setup
dump=DumpForPlots
base=${ben_arch}_${sample}

####################
## Run Benchmarks ##
####################

## compile with appropriate options
make distclean ${mOpt}
make ${mOpt}

## Vectorization Benchmarks
for nvu in "${nvus[@]}"
do
    make clean ${mOpt}
    make ${mOpt} USE_INTRINSICS:=-DMPT_SIZE=${nvu}

    ## Common base executable
    oBase=${base}
    bExe="${exe} --fit-std-only --num-thr ${minth} --num-events ${nevents}"

    ## Building-only benchmark
    echo "${oBase}: Benchmark [nTH:${minth}, nVU:${nvu}]"
    ${bExe} >& log_${oBase}_NVU${nvu}_NTH${minth}.txt
done

## Final cleanup
make distclean ${mOpt}
