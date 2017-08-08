#!/bin/bash

dir=${1:-cmsswval-pu70}
base=SNB_CMSSW_PU70

echo "Moving plots and text files locally to ${dir}"
for seed in SimSeed CMSSeed
do
    outdir=${dir}/${seed}
    mkdir -p ${outdir}
 
    mv ${base}_${seed}_*.png ${outdir}
    for build in BH STD CE
    do
	vbase=validation_${base}_${seed}_${build}
	mv ${vbase}/totals_${vbase}.txt ${outdir}
    done
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
for seed in SimSeed CMSSeed
do
    for build in BH STD CE
    do
	testbase=${base}_${seed}_${build}
	rm -rf validation_${testbase}
	rm -rf log_${testbase}_NVU8int_NTH24_val.txt 
    done
done

rm -rf ${dir}
