#! /bin/bash

## In the case this is run separately from main script
[ -z "$ROOTSYS" ] && source ~matevz/root/bin/thisroot.sh
source xeon_scripts/common_variables.sh

## Common setup
dir=/data/nfsmic/slava77/samples/2017/pass-4874f28/initialStep
subdir=PU70HS/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU/a
file=memoryFile.fv3.clean.writeAll.recT.072617.bin
nevents=500

maxth=24
maxvu=8
exe="./mkFit/mkFit --cmssw-n2seeds --num-thr ${maxth} --input-file ${dir}/${subdir}/${file} --num-events ${nevents}"

## validation function
doVal()
{
    bN=${1}
    bO=${2}
    vN=${3}
    vO=${4}

    oBase=${val_arch}_${sample}_${bN}
    bExe="${exe} ${vO} --build-${bO}"
    
    echo "${oBase}: ${vN} [nTH:${maxth}, nVU:${maxvu}int]"
    ${bExe} >& log_${oBase}_NVU${maxvu}int_NTH${maxth}_${vN}.txt
    
    mv valtree.root valtree_${oBase}_${vN}.root
}		

## Compile once
make clean
mVal="WITH_ROOT:=yes"
make -j 12 ${mVal}

## Special simtrack validation vs cmssw tracks
doVal CMSSW cmssw SIMVAL "--sim-val-for-cmssw --read-cmssw-tracks"

## SIMVAL == sim tracks as reference, CMSSWVAL == cmssw tracks as reference
for vV in "SIMVAL --sim-val" "CMSSWVAL --cmssw-val-hit"
do echo ${vV} | while read -r vN vO
    do
	for bV in "BH bh" "STD std" "CE ce" "FV fv"
	do echo ${bV} | while read -r bN bO
	    do
		doVal ${bN} ${bO} ${vN} "${vO} --backward-fit"
	    done
	done
    done
done

## clean up
make clean ${mVal}

## Compute observables and make images
for vV in "SIMVAL 0" "CMSSWVAL 1"
do echo ${vV} | while read -r vN vO
    do
	tbase=${val_arch}_${sample}
	for build in BH STD CE FV CMSSW
	do
	    if [[ "${build}" == "CMSSW" ]] && [[ "${vN}" == "CMSSWVAL" ]] ; then
		continue
	    fi

	    echo "Computing observables for: ${tbase} ${build} ${vN}"
	    root -b -q -l plotting/runValidation.C\(\"_${tbase}_${build}_${vN}\",0,${vO}\)
	done

	## overlay histograms
	echo "Overlaying histograms for: ${tbase} ${vN}"
	root -b -q -l plotting/makeValidation.C\(\"${tbase}\",\"_${vN}\",${vO}\)
    done
done

make distclean ${mVal}
