#! /bin/bash

# remove tmp dir on KNL remotely
echo "Removing tmp dir on KNL remotely"
ssh ${KNL_HOST} bash -c "'
rm -rf ${KNL_WORKDIR}/${KNL_TEMPDIR}
exit
'"
