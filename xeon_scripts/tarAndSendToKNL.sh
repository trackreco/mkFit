#! /bin/bash

# in case this is sent separately
source xeon_scripts/common_variables.sh

# tar up the directory
echo "Tarring directory for KNL... make sure it is clean!"
repo=mictest.tar.gz
tar --exclude-vcs --exclude='*.gz' --exclude='validation*' --exclude='*.root' --exclude='log_*' --exclude='*.png' --exclude='*.o' --exclude='*.om' --exclude='*.d' --exclude='*.optrpt' -zcvf  ${repo} *

# mkdir tmp dir on KNL remotely
echo "Making tmp dir on KNL remotely"
SSHO ${KNL_HOST} bash -c "'
mkdir -p ${KNL_WORKDIR}/${KNL_TEMPDIR}
exit
'"

# copy tarball
echo "Copying tarball to KNL"
scp ${repo} ${KNL_HOST}:${KNL_WORKDIR}/${KNL_TEMPDIR}

# unzip tarball remotely
echo "Untarring repo on KNL remotely"
SSHO ${KNL_HOST} bash -c "'
cd ${KNL_WORKDIR}/${KNL_TEMPDIR}
tar -zxvf ${repo}
rm ${repo}
'"

# remove local tarball
echo "Remove local repo tarball"
rm ${repo}
