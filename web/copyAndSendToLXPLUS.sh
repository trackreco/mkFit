#! /bin/bash

# in case this is run alone
source xeon_scripts/common_variables.sh

dir=${1:-benchmarks}

# first copy all the plots
scp -r ${dir} ${LXPLUS_HOST}:${LXPLUS_WORKDIR}/${LXPLUS_OUTDIR}

# Make outdir nice and pretty
echo "Executing remotely ./makereadable.sh ${dir}"
ssh ${LXPLUS_HOST} bash -c "'
cd ${LXPLUS_OUTDIR}
./makereadable.sh ${dir}
exit
'"
