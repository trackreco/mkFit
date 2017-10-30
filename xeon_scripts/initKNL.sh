#! /bin/bash

echo "Initializing KNL setup"
source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/git/1.8.3.1/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/gcc/6.3.0/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/lcg/root/6.08.07/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/pcre/8.37/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/xz/5.2.2/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/libtiff/4.0.3/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/libpng/1.6.16/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/libjpeg-turbo/1.3.1/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/sqlite/3.15.1/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/gdb/7.12/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/python/2.7.11/etc/profile.d/init.sh

. /opt/intel/bin/compilervars.sh -arch intel64
export TBB_PREFIX=$TBBROOT
