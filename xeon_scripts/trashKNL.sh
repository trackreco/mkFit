#! /bin/bash

# in case this is called separately
source xeon_scripts/common_variables.sh

# remove tmp dir on KNL remotely
echo "Removing tmp dir on KNL remotely"
SSHO ${KNL_HOST} bash -c "'
rm -rf ${KNL_WORKDIR}/${KNL_TEMPDIR}
exit
'"
