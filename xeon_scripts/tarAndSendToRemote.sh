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
else 
    echo ${remote_arch} "is not a valid architecture! Exiting..."
    exit
fi

##################
## Tar and Send ##
##################

# tar up the directory
echo "Tarring directory for ${remote_arch}... make sure it is clean!"
repo=mictest.tar.gz
tar --exclude-vcs --exclude='*.gz' --exclude='validation*' --exclude='*.root' --exclude='log_*' --exclude='*.png' --exclude='*.o' --exclude='*.om' --exclude='*.d' --exclude='*.optrpt' -zcvf  ${repo} *

# mkdir tmp dir on remote arch
echo "Making tmp dir on ${remote_arch} remotely"
SSHO ${HOST} bash -c "'
mkdir -p ${DIR}
exit
'"

# copy tarball
echo "Copying tarball to ${remote_arch}"
scp ${repo} ${HOST}:${DIR}

# unzip tarball remotely
echo "Untarring repo on ${remote_arch} remotely"
SSHO ${HOST} bash -c "'
cd ${DIR}
tar -zxvf ${repo}
rm ${repo}
'"

# remove local tarball
echo "Remove local repo tarball"
rm ${repo}
