#! /bin/bash

###########
## Input ##
###########

suite=${1:-"forPR"} # which set of benchmarks to run: full, forPR, forConf

###################
## Configuration ##
###################

source xeon_scripts/common-variables.sh ${suite}
source xeon_scripts/init-env.sh

## Common file setup
dir=/data2/slava77/samples/2017/pass-c93773a/initialStep
subdir=PU70HS/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU
file=memoryFile.fv3.clean.writeAll.recT.082418-25daeda.bin
nevents=500

## Common executable setup
maxth=64
maxvu=16
maxev=32
seeds="--cmssw-n2seeds"
exe="./mkFit/mkFit --silent ${seeds} --num-thr ${maxth} --num-thr-ev ${maxev} --input-file ${dir}/${subdir}/${file} --num-events ${nevents}"

## Common output setup
tmpdir="tmp"
base=${val_arch}_${sample}

## flag to save sim info for matched tracks since track states not read in
siminfo="--try-to-save-sim-info"

## backward fit flag
bkfit="--backward-fit-pca"

## validation options: SIMVAL == sim tracks as reference, CMSSWVAL == cmssw tracks as reference
SIMVAL="SIMVAL --sim-val ${siminfo} ${bkfit}"
CMSSWVAL="CMSSWVAL --cmssw-val-fhit-bprm ${bkfit}"
declare -a vals=(SIMVAL CMSSWVAL)

## plotting options
SIMPLOT="SIMVAL 0"
CMSSWPLOT="CMSSWVAL 1"
declare -a plots=(SIMPLOT CMSSWPLOT)

## special cmssw dummy build
CMSSW="CMSSW cmssw SIMVAL --sim-val-for-cmssw ${siminfo} --read-cmssw-tracks"

###############
## Functions ##
###############

## validation function
function doVal()
{
    local bN=${1}
    local bO=${2}
    local vN=${3}
    local vO=${4}

    local oBase=${val_arch}_${sample}_${bN}
    local bExe="${exe} ${vO} --build-${bO}"
    
    echo "${oBase}: ${vN} [nTH:${maxth}, nVU:${maxvu}int, nEV:${maxev}]"
    ${bExe} >& log_${oBase}_NVU${maxvu}int_NTH${maxth}_NEV${maxev}_${vN}.txt
    
    # hadd output files for this test, then move to temporary directory
    hadd -O valtree.root valtree_*.root
    rm valtree_*.root
    mv valtree.root ${tmpdir}/valtree_${oBase}_${vN}.root
}		

## plotting function
function plotVal()
{
    local base=${1}
    local bN=${2}
    local pN=${3}
    local pO=${4}

    echo "Computing observables for: ${base} ${bN} ${pN}"
    root -b -q -l plotting/runValidation.C\(\"_${base}_${bN}_${pN}\",${pO}\)
}

########################
## Run the validation ##
########################

## Compile once
make clean
mVal="-j 32 WITH_ROOT:=1 AVX_512:=1"
make ${mVal}
mkdir -p ${tmpdir}

## Special simtrack validation vs cmssw tracks
echo ${CMSSW} | while read -r bN bO vN vO
do
    doVal "${bN}" "${bO}" "${vN}" "${vO}"
done

## Run validation for standard build options
for val in "${vals[@]}"
do echo ${!val} | while read -r vN vO
    do
	for build in "${val_builds[@]}"
	do echo ${!build} | while read -r bN bO
	    do
		doVal "${bN}" "${bO}" "${vN}" "${vO}"
	    done
	done
    done
done

## clean up
make clean ${mVal}
mv tmp/valtree_*.root .
rm -rf ${tmpdir}

## Compute observables and make images
for plot in "${plots[@]}"
do echo ${!plot} | while read -r pN pO
    do
        ## Compute observables for special dummy CMSSW
	if [[ "${pN}" == "SIMVAL" ]]
	then
	    echo ${CMSSW} | while read -r bN bO val_extras
	    do
		plotVal "${base}" "${bN}" "${pN}" "${pO}"
	    done
	fi

	## Compute observables for builds chosen 
	for build in "${val_builds[@]}"
	do echo ${!build} | while read -r bN bO
	    do
		plotVal "${base}" "${bN}" "${pN}" "${pO}"
	    done
	done
	
	## overlay histograms
	echo "Overlaying histograms for: ${base} ${vN}"
	root -b -q -l plotting/makeValidation.C\(\"${base}\",\"_${pN}\",${pO},\"${suite}\"\)
    done
done

make distclean ${mVal}
