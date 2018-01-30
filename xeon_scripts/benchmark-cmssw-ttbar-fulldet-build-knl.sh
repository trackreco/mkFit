#! /bin/bash

# in case this is called separately
source xeon_scripts/common_variables.sh

# execute knl tests remotely
echo "Executing KNL tests remotely..."
SSHO ${KNL_HOST} bash -c "'
cd ${KNL_WORKDIR}/${KNL_TEMPDIR}
source xeon_scripts/initKNL.sh
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build.sh KNL
exit
'"

# copy logs back for plotting
echo "Copying logs back from KNL for plotting"
scp ${KNL_HOST}:${KNL_WORKDIR}/${KNL_TEMPDIR}/log_KNL_${physics_sample}_*.txt .

# destroy tmp files
./xeon_scripts/trashKNL.sh
