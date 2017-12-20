#! /bin/bash

# remove files remotely
echo "Removing files on KNC remotely"
ssh < /dev/null -o 'StrictHostKeyChecking no'  ${KNC_HOST} bash -c "'
cd ${KNC_WORKDIR}
rm -rf Geoms/
rm -rf mkFit/
rm -rf *.so
exit
'"
