#! /bin/bash

########################
## Command Line Input ##
########################

remote_arch=${1} # SNB, KNL, SKL-SP
suite=${2:-"forPR"} # which set of benchmarks to run: full, forPR, forConf

###################
## Configuration ##
###################

source xeon_scripts/common-variables.sh ${suite}
source xeon_scripts/init-env.sh

# architecture dependent settings
if [[ "${remote_arch}" == "SNB" ]]
then
    HOST=${SNB_HOST}
    DIR=${SNB_WORKDIR}/${SNB_TEMPDIR}
elif [[ "${remote_arch}" == "KNL" ]]
then
    HOST=${KNL_HOST}
    DIR=${KNL_WORKDIR}/${KNL_TEMPDIR}
elif [[ "${remote_arch}" == "LNX" ]]
then 
    HOST=${LNX_HOST}
    DIR=${LNX_WORKDIR}/${LNX_TEMPDIR}
elif [[ "${remote_arch}" == "LNX-S" ]]
then
    HOST=${LNXS_HOST}
    DIR=${LNXS_WORKDIR}/${LNXS_TEMPDIR}
else 
    echo ${remote_arch} "is not a valid architecture! Exiting..."
    exit
fi

###################
## Run The Tests ##
###################

# execute tests remotely
echo "Executing ${remote_arch} tests remotely..."
SSHO ${HOST} bash -c "'
echo "Executing test ${remote_arch} tests remotely..."
cd ${DIR}
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build.sh ${remote_arch} ${suite}
exit
'"

# copy logs back for plotting
echo "Copying logs back from ${remote_arch} for plotting"
scp ${HOST}:${DIR}/log_${remote_arch}_${sample}_*.txt .

# destroy tmp files
echo "Removing tmp dir on ${remote_arch} remotely"
SSHO ${HOST} bash -c "'
rm -rf ${DIR}
exit
'"
