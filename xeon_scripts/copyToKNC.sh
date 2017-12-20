#! /bin/bash

# make temp dirs
echo "Making tmp dirs on KNC remotely"
ssh < /dev/null -o StrictHostKeyChecking=no ${KNC_HOST} bash -c "'
cd ${KNC_WORKDIR}
mkdir -p Geoms
mkdir -p mkFit
exit
'"

# copy the executables + shared object
scp libMicCore-mic.so ${KNC_HOST}:${KNC_WORKDIR}
scp mkFit/mkFit-mic ${KNC_HOST}:${KNC_WORKDIR}/mkFit
scp Geoms/CMS-2017-mic.so ${KNC_HOST}:${KNC_WORKDIR}/Geoms
scp Geoms/CylCowWLids-mic.so ${KNC_HOST}:${KNC_WORKDIR}/Geoms
