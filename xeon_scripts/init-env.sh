#!/bin/bash

source /cvmfs/cms.cern.ch/slc7_amd64_gcc820/lcg/root/6.18.04-bcolbf/etc/profile.d/init.sh
# workaround for https://github.com/cms-sw/cmsdist/issues/5574
# remove when we switch to a ROOT build where that issues is fixed
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBJPEG_TURBO_ROOT/lib64
source /opt/intel/bin/compilervars.sh intel64
