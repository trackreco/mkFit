#!/bin/bash

dir=${1:-toymcval}
base=SNB_ToyMC_FullDet

echo "Moving plots and text files locally to ${dir}"
mkdir -p ${dir}
mv ${base}_*.png ${dir}
for build in BH STD CE
do
    vbase=validation_${base}_${build}
    mv ${vbase}/totals_${vbase}.txt ${dir}
done

host=kmcdermo@lxplus.cern.ch
whost=${host}":~/www"
echo "Moving plots and text files remotly to ${whost}"
scp -r ${dir} ${whost}

echo "Executing remotely ./makereadable.sh ${dir}"
ssh ${host} bash -c "'
cd www
./makereadable.sh ${dir}
exit
'"

echo "Removing local files"
for build in BH STD CE
do
    testbase=${base}_${build}
    rm -rf validation_${testbase}
    rm -rf log_${testbase}_NVU8int_NTH24_val.txt 
done

rm -rf ${dir}
