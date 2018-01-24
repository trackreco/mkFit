#! /bin/bash

## SNB or KNL
arch=${1}

## In the case this is run separately from main script
source xeon_scripts/common_variables.sh

## Common setup
subdir=2017/pass-4874f28/initialStep/PU70HS/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU/a/
file=memoryFile.fv3.clean.writeAll.recT.072617.bin
nevents=20
minth=1
minvu=1
dump=DumpForPlots
seeds="--cmssw-n2seeds"

## Platform specific settings
if [ ${arch} == "SNB" ]
then
    mOpt="-j 12"
    dir=/data/nfsmic/slava77/samples
    base=${arch}_${physics_sample}
    maxth=24
    maxvu=8
    exe="./mkFit/mkFit ${seeds} --input-file ${dir}/${subdir}/${file}"
    declare -a nths=("1" "2" "4" "6" "8" "12" "16" "20" "24")
    declare -a nvus=("1" "2" "4" "8")
    declare -a nevs=("1" "2" "4" "8" "12")
elif [ ${arch} == "KNL" ]
then
    mOpt="-j 64 AVX_512:=1"
    dir=/data1/work/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000
    base=${arch}_${physics_sample}
    maxth=256
    maxvu=16
    exe="./mkFit/mkFit ${seeds} --input-file ${dir}/${subdir}/${file}"
    declare -a nths=("1" "2" "4" "8" "16" "32" "64" "96" "128" "160" "192" "224" "256")
    declare -a nvus=("1" "2" "4" "8" "16")
    declare -a nevs=("1" "2" "4" "8" "16" "32" "64" "128")
else 
    echo ${arch} "is not a valid architecture! Exiting..."
    exit
fi

## compile with appropriate options
make distclean ${mOpt}
make ${mOpt}

## Parallelization Benchmarks
for nth in "${nths[@]}"
do
    for bV in "BH bh" "STD std" "CE ce" "FV fv"
    do echo ${bV} | while read -r bN bO
	do
	    ## Base executable
	    oBase=${base}_${bN}
	    bExe="${exe} --build-${bO} --num-thr ${nth}"

	    ## Building-only benchmark
	    echo "${oBase}: Benchmark [nTH:${nth}, nVU:${maxvu}int]"
	    ${bExe} --num-events ${nevents} >& log_${oBase}_NVU${maxvu}int_NTH${nth}.txt

	    ## Multiple Events in Flight benchmark
	    if [ "${bN}" == "CE" ] || [ "${bN}" == "FV" ]
	    then
		for nev in "${nevs[@]}"
		do
		    if (( ${nev} <= ${nth} ))
		    then
			nproc=$(( ${nevents} * ${nev} ))
			echo "${oBase}: Benchmark [nTH:${nth}, nVU:${maxvu}int, nEV:${nev}]"
			${bExe} --silent --num-thr-ev ${nev} --num-events ${nproc} >& log_${oBase}_NVU${maxvu}int_NTH${nth}_NEV${nev}.txt
		    fi
		done
	    fi

	    ## nHits validation
	    if (( ${nth} == ${maxth} ))
	    then
		echo "${oBase}: Text dump for plots [nTH:${nth}, nVU:${maxvu}int]"
		${bExe} --dump-for-plots --quality-val --read-cmssw-tracks --num-events ${nevents} >& log_${oBase}_NVU${maxvu}int_NTH${nth}_${dump}.txt
	    fi
	done
    done
done

## Vectorization Benchmarks
for nvu in "${nvus[@]}"
do
    make clean ${mOpt}
    make ${mOpt} USE_INTRINSICS:=-DMPT_SIZE=${nvu}

    for bV in "BH bh" "STD std" "CE ce" # "FV fv"
    do echo ${bV} | while read -r bN bO
	do
	    ## Common base executable
	    oBase=${base}_${bN}
	    bExe="${exe} --build-${bO} --num-thr ${minth} --num-events ${nevents}"

	    ## Building-only benchmark
	    echo "${oBase}: Benchmark [nTH:${minth}, nVU:${nvu}]"
	    ${bExe} >& log_${oBase}_NVU${nvu}_NTH${minth}.txt

	    ## nHits validation
	    if (( ${nvu} == ${minvu} ))
	    then
		echo "${oBase}: Text dump for plots [nTH:${minth}, nVU:${nvu}]"
		${bExe} --dump-for-plots --quality-val --read-cmssw-tracks >& log_${oBase}_NVU${nvu}_NTH${minth}_${dump}.txt
	    fi
	done
    done
done

## Final cleanup
make distclean ${mOpt}
