#! /bin/bash

## In the case this is run separately from main script
[ -z "$ROOTSYS" ] && source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/lcg/root/6.12.07-gnimlf/etc/profile.d/init.sh
source xeon_scripts/common_variables.sh

## Common setup
dir=/data2/slava77/samples/2017/pass-4874f28/initialStep
subdir=PU70HS/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU/a
file=memoryFile.fv3.clean.writeAll.recT.072617.bin
nevents=500
tmpdir="tmp"

maxth=64
maxvu=16
#maxev=12 FIXME: comment out MEIF stuff for now
maxev=1
exe="./mkFit/mkFit --cmssw-n2seeds --num-thr ${maxth} --num-thr-ev ${maxev} --input-file ${dir}/${subdir}/${file} --num-events ${nevents}"

## validation function
doVal()
{
    bN=${1}
    bO=${2}
    vN=${3}
    vO=${4}

    oBase=${val_arch}_${sample}_${bN}
    bExe="${exe} ${vO} --build-${bO}"
    
    echo "${oBase}: ${vN} [nTH:${maxth}, nVU:${maxvu}int, nEV:${maxev}]"
    ${bExe} >& log_${oBase}_NVU${maxvu}int_NTH${maxth}_NEV${maxev}_${vN}.txt
    
#    FIXME: comment out MEIF stuff for now
#    hadd valtree.root valtree_*.root
#    rm valtree_*.root
    mv valtree.root ${tmpdir}/valtree_${oBase}_${vN}.root
}		

## Compile once
make clean
mVal="WITH_ROOT:=yes AVX_512:=1"
make -j 32 ${mVal}
mkdir -p ${tmpdir}

## Special simtrack validation vs cmssw tracks
doVal CMSSW cmssw SIMVAL "--sim-val-for-cmssw --read-cmssw-tracks"

## SIMVAL == sim tracks as reference, CMSSWVAL == cmssw tracks as reference
for vV in "SIMVAL --sim-val" "CMSSWVAL --cmssw-val-fhit-bprm"
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
mv tmp/valtree_*.root .
rm -rf ${tmpdir}

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
	    root -b -q -l plotting/runValidation.C\(\"_${tbase}_${build}_${vN}\",${vO}\)
	done

	## overlay histograms
	echo "Overlaying histograms for: ${tbase} ${vN}"
	root -b -q -l plotting/makeValidation.C\(\"${tbase}\",\"_${vN}\",${vO}\)
    done
done

make distclean ${mVal}
