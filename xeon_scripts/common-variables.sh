#! /bin/bash

# command line input
suite=${1:-"forPR"} # which set of benchmarks to run: full, forPR, forConf

# samples
export sample=CMSSW_TTbar_PU70

# Validation architecture
export val_arch=SKL-SP

# vars for KNL
export KNL_HOST=${USER}@phi2.t2.ucsd.edu
export KNL_WORKDIR=/data1/work/${USER}
export KNL_TEMPDIR=tmp

# vars for SNB
export SNB_HOST=${USER}@phi1.t2.ucsd.edu
export SNB_WORKDIR=/data2/nfsmic/${USER}
export SNB_TEMPDIR=tmp

# vars for LXPLUS
export LXPLUS_HOST=${USER}@lxplus.cern.ch
export LXPLUS_WORKDIR=/afs/cern.ch/user/${USER:0:1}/${USER}
export LXPLUS_OUTDIR=www

# SSH options
function SSHO()
{
    ssh -o StrictHostKeyChecking=no < /dev/null "$@"
}
export -f SSHO

#################
## Build Types ##
#################

export BH="BH bh"
export STD="STD std"
export CE="CE ce"
export FV="FV fv"

# which set of builds to use based on input from command line
if [[ "${suite}" == "full" ]]
then
    declare -a ben_builds=(BH STD CE FV)
    declare -a val_builds=(BH STD CE FV)
elif [[ "${suite}" == "forPR" ]]
then
    declare -a ben_builds=(BH CE)
    declare -a val_builds=(STD CE)
elif [[ "${suite}" == "forConf" ]]
then
    declare -a ben_builds=(CE)
    declare -a val_builds=(CE)
else
    echo ${suite} "is not a valid benchmarking suite option! Exiting..."
    exit
fi

# set dependent arrays
th_builds=() ## for parallelization tests
vu_builds=() ## for vectorization tests
meif_builds=() ## for multiple-events-in-flight tests
text_builds=() ## for text dump comparison tests

# loop over ben_builds and set dependent arrays, export when done
for build in "${ben_builds[@]}"
do
    # set th builds : all benchmarks!
    th_builds+=("${build}")
    
    # set vu builds : exclude FV since it does not have a meaningful implementation outside of max VU
    if [[ "${build}" != "FV" ]]
    then
	vu_builds+=("${build}")
    fi
    
    # set meif builds : only do CE and FV
    if [[ "${build}" == "CE" ]] || [[ "${build}" == "FV" ]]
    then
	meif_builds+=("${build}")
    fi
done

# set text dump builds: need builds matched in both TH and VU tests
for th_build in "${th_builds[@]}"
do 
    for vu_build in "${vu_builds[@]}"
    do
	if [[ "${th_build}" == "${vu_build}" ]]
	then
	    text_builds+=("${th_build}")
	fi
    done
done
for vu_build in "${vu_builds[@]}"
do 
    for th_build in "${th_builds[@]}"
    do
	if [[ "${vu_build}" == "${th_build}" ]]
	then
	    text_builds+=("${vu_build}")
	fi
    done
done    
text_builds=( $( for build in "${text_builds[@]}"; do echo "${build}"; done | sort | uniq | xargs ) )

export ben_builds val_builds th_builds vu_builds meif_builds text_builds

# meif checking
function CheckIfMEIF ()
{
    local build=${1}
    local result="false"

    for meif_build in "${meif_builds[@]}"
    do 
	if [[ "${meif_build}" == "${build}" ]]
	then
	    result="true"
	    break
	fi
    done
    
    echo "${result}"
}
export -f CheckIfMEIF

# text checking
function CheckIfText ()
{
    local build=${1}
    local result="false"

    for text_build in "${text_builds[@]}"
    do 
	if [[ "${text_build}" == "${build}" ]]
	then
	    result="true"
	    break
	fi
    done

    echo "${result}"
}
export -f CheckIfText
