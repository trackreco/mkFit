#! /bin/bash

# samples
export rootval_sample=CMSSW_10mu
export physics_sample=CMSSW_TTbar_PU70

# vars for KNL
export KNL_HOST=${USER}@phi2.t2.ucsd.edu
export KNL_WORKDIR=/data1/work/${USER}
export KNL_TEMPDIR=tmp

# vars for LXPLUS
export LXPLUS_HOST=${USER}@lxplus.cern.ch
export LXPLUS_WORKDIR=/afs/cern.ch/user/${USER:0:1}/${USER}
export LXPLUS_OUTDIR=www

# SSH options
SSHO() {
    ssh -o StrictHostKeyChecking=no < /dev/null "$@"
}
export -f SSHO

# Validation architecture
export val_arch=SNB
