#! /bin/bash

# in case this is sent separately
source xeon_scripts/common_variables.sh

# tar up the directory
echo "Tarring directory for SNB... make sure it is clean!"
repo=mictest.tar.gz
tar --exclude-vcs --exclude='*.gz' --exclude='validation*' --exclude='*.root' --exclude='log_*' --exclude='*.png' --exclude='*.o' --exclude='*.om' --exclude='*.d' --exclude='*.optrpt' -zcvf  ${repo} *

# mkdir tmp dir on SNB remotely
echo "Making tmp dir on SNB remotely"
SSHO ${SNB_HOST} bash -c "'
mkdir -p ${SNB_WORKDIR}/${SNB_TEMPDIR}
exit
'"

# copy tarball
echo "Copying tarball to SNB"
scp ${repo} ${SNB_HOST}:${SNB_WORKDIR}/${SNB_TEMPDIR}

# unzip tarball remotely
echo "Untarring repo on SNB remotely"
SSHO ${SNB_HOST} bash -c "'
cd ${SNB_WORKDIR}/${SNB_TEMPDIR}
tar -zxvf ${repo}
rm ${repo}
'"

# remove local tarball
echo "Remove local repo tarball"
rm ${repo}
