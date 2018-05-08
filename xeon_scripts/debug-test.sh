#! /bin/bash

## Run this script on SNB before benchmarking and submitting a PR to ensure your new commit does not break the debug mode!

## In the case this is run separately from main script
[ -z "$ROOTSYS" ] && source ~matevz/root/bin/thisroot.sh
source xeon_scripts/common_variables.sh

## Common setup
dir=/data/nfsmic/slava77/samples/2017/pass-4874f28/initialStep
subdir=10muPt0p5to10HS
file=memoryFile.fv3.clean.writeAll.recT.072617.bin

## config for debug
nevents=10
maxth=1
maxvu=1
maxev=1

## base executable
exe="./mkFit/mkFit --cmssw-n2seeds --num-thr ${maxth} --num-thr-ev ${maxev} --input-file ${dir}/${subdir}/${file} --num-events ${nevents}"

## Compile once
mOpt="DEBUG:=yes WITH_ROOT:=yes USE_INTRINSICS:=-DMPT_SIZE=${maxvu}"
make distclean ${mOpt}
make -j 12 ${mOpt}

## test each build routine to be sure it works!
for bV in "BH bh" "STD std" "CE ce" "FV fv"
do echo ${bV} | while read -r bN bO
    do
	oBase=${val_arch}_${sample}_${bN}
	bExe="${exe} --build-${bO}"

	echo "${oBase}: ${vN} [nTH:${maxth}, nVU:${maxvu}, nEV:${maxev}]"
	${bExe} >& log_${oBase}_NVU${maxvu}_NTH${maxth}_NEV${maxev}_"DEBUG".txt
    done
done

## clean up
make distclean ${mOpt}
