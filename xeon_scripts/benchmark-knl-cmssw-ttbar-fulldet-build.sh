#! /bin/bash

# vars for directories
host=${USER}@phi2.t2.ucsd.edu
workdir=/data1/work/${USER}
tmpdir=tmp

# execute knl tests remotely
echo "Executing KNL tests remotely..."
ssh ${host} bash -c "'
cd ${workdir}/${tmpdir}
source xeon_scripts/initKNL.sh
./xeon_scripts/benchmark-knl-cmssw-ttbar-fulldet-build-remote.sh
exit
'"

# copy logs back for plotting
echo "Copying logs back from KNL for plotting"
scp ${host}:${workdir}/${tmpdir}/log_*.txt .
