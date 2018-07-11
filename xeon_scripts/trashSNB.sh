#! /bin/bash

# in case this is called separately
source xeon_scripts/common-variables.sh
source xeon_scripts/init-env.sh

# remove tmp dir on SNB remotely
echo "Removing tmp dir on SNB remotely"
SSHO ${SNB_HOST} bash -c "'
rm -rf ${SNB_WORKDIR}/${SNB_TEMPDIR}
exit
'"
