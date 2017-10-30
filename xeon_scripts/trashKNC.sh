#! /bin/bash

# vars for directories
host=${USER}@phiphi-mic0.t2.ucsd.edu
workdir=/home/${USER}

# remove files remotely
echo "Removing files on KNC remotely"
ssh ${host} bash -c "'
cd ${workdir}
rm -rf Geoms/
rm -rf mkFit/
rm -rf *.so
exit
'"
