#!/bin/bash

dir=${1:-cmsswval}
base=SNB_CMSSW_10mu

echo "Moving plots and text files locally to ${dir}"
for seed in SimSeed CMSSeed
do
    for region in ECN2 ECN1 BRL ECP1 ECP2 FullDet
    do
	outdir=${dir}/${seed}/${region}
	mkdir -p ${outdir}
 
	srbase=${seed}_${region}
	mv ${base}_${srbase}_*.png ${outdir}
	for build in BH STD CE
	do
	    vbase=validation_${base}_${srbase}_${build}
	    mv ${vbase}/totals_${vbase}.txt ${outdir}
	done
    done
    sdir=${dir}/${seed}
    mv ${sdir}/FullDet/*png ${sdir}/FullDet/*txt ${sdir}
    rm -rf ${sdir}/FullDet
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
    for region in ECN2 ECN1 BRL ECP1 ECP2 FullDet
    do
	srbase=${seed}_${region}
	for build in BH STD CE
	do
	    testbase=${base}_${srbase}_${build}
	    rm -rf validation_${testbase}
	    rm -rf log_${testbase}_NVU8int_NTH24_val.txt 
	done
    done
done

rm -rf ${dir}
