#! /bin/bash

# tar up the directory
echo "Tarring directory for KNL... make sure it is clean!"
repo=mictest.tar.gz
tar -zcvf ${repo} *

# vars for directories
host=${USER}@phi2.t2.ucsd.edu
workdir=/data1/work/${USER}
tmpdir=tmp

# mkdir tmp dir on KNL remotely
echo "Making tmp dir on KNL remotely"
ssh ${host} bash -c "'
cd ${workdir}
mkdir -p ${tmpdir}
exit
'"

# copy tarball
echo "Copying tarball to KNL"
scp ${repo} ${host}:${workdir}/${tmpdir}

# unzip tarball remotely
echo "Untarring repo on KNL remotely"
ssh ${host} bash -c "'
cd ${workdir}/${tmpdir}
tar -zxvf ${repo}
rm ${repo}
'"

echo "Remove local repo tarball"
rm ${repo}
