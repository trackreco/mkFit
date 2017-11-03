#! /bin/bash

# sample
export sample=CMSSW_TTbar_PU70

# vars for KNC
export KNC_HOST=${USER}@phiphi-mic0.t2.ucsd.edu
export KNC_WORKDIR=/home/${USER}

# vars for KNL
export KNL_HOST=${USER}@phi2.t2.ucsd.edu
export KNL_WORKDIR=/data1/work/${USER}
export KNL_TEMPDIR=tmp

# vars for LXPLUS
export LXPLUS_HOST=${USER}@lxplus.cern.ch
export LXPLUS_WORKDIR=/afs/cern.ch/user/${USER:0:1}/${USER}
export LXPLUS_OUTDIR=www