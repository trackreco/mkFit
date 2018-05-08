#! /bin/bash

echo "Initializing SNB setup"
source /opt/rh/devtoolset-2/enable;
source /opt/intel/bin/iccvars.sh intel64;
source /opt/intel/vtune_amplifier_xe/amplxe-vars.sh
source /opt/intel/bin/compilervars.sh -arch intel64
