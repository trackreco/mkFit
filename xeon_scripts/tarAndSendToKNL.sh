#! /bin/bash

# tar up the directory
echo "Tarring directory for KNL... make sure it is clean!"
repo=mictest.tar.gz
tar -zcvf ${repo} *

# mkdir tmp dir on KNL remotely
echo "Making tmp dir on KNL remotely"
ssh ${KNL_HOST} bash -c "'
mkdir -p ${KNL_WORKDIR}/${KNL_TEMPDIR}
exit
'"

# copy tarball
echo "Copying tarball to KNL"
scp ${repo} ${KNL_HOST}:${KNL_WORKDIR}/${KNL_TEMPDIR}

# unzip tarball remotely
echo "Untarring repo on KNL remotely"
ssh ${KNL_HOST} bash -c "'
cd ${KNL_WORKDIR}/${KNL_TEMPDIR}
tar -zxvf ${repo}
rm ${repo}
'"

# remove local tarball
echo "Remove local repo tarball"
rm ${repo}
