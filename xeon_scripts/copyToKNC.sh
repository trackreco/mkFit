#! /bin/bash

# vars for directories
host=${USER}@phiphi-mic0.t2.ucsd.edu
workdir=/home/${USER}

# make temp dirs
echo "Making tmp dirs on KNC remotely"
ssh ${host} bash -c "'
cd ${workdir}
mkdir -p Geoms
mkdir -p mkFit
exit
'"

# copy the executables + shared object
scp libMicCore-mic.so ${host}:${workdir}
scp mkFit/mkFit-mic ${host}:${workdir}/mkFit
scp Geoms/CMS-2017-mic.so ${host}:${workdir}/Geoms
scp Geoms/CylCowWLids-mic.so ${host}:${workdir}/Geoms
