#! /bin/bash

# vars for directories
host=${USER}@phi2.t2.ucsd.edu
workdir=/data1/work/${USER}
tmpdir=tmp

# remove tmp dir on KNL remotely
echo "Removing tmp dir on KNL remotely"
ssh ${host} bash -c "'
rm -rf ${workdir}/${tmpdir}
exit
'"
