#! /bin/bash

# command line input
dir=${1:-"benchmarks"} # Main output dir name
suite=${2:-"forPR"} # which set of benchmarks to run: full, forPR, forConf

# in case this is run alone
source xeon_scripts/common-variables.sh ${suite}
source xeon_scripts/init-env.sh

# first tar the directory to be sent
echo "Tarring plot directory"
tarball=${dir}.tar.gz
tar -zcvf ${tarball} ${dir}

# then send it!
scp -r ${tarball} ${LXPLUS_HOST}:${LXPLUS_WORKDIR}/${LXPLUS_OUTDIR}

# Make outdir nice and pretty
echo "Unpack tarball and execute remotely ./makereadable.sh ${dir}"
SSHO ${LXPLUS_HOST} bash -c "'
cd ${LXPLUS_OUTDIR}
tar -zxvf ${tarball}
./makereadable.sh ${dir}
rm -rf ${tarball}
exit
'"

# remove local tarball
echo "Remove local tarball of plots"
rm ${tarball}
