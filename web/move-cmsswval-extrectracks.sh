#!/bin/bash

dir=${1:-cmsswval-ext}
base=SNB_CMSSW_10mu

echo "Moving plots and text files locally to ${dir}"
for region in ECN2 ECN1 BRL ECP1 ECP2 FullDet
do
    outdir=${dir}/${region}
    mkdir -p ${outdir}
    
    mv ${base}_${region}_*.png ${outdir}
    for build in BH STD CE
    do
	vbase=validation_${base}_${region}_${build}
	mv ${vbase}/totals_${vbase}_cmssw.txt ${outdir}
    done
done

mv ${dir}/FullDet/*png ${dir}/FullDet/*txt ${dir}
rm -rf ${dir}/FullDet

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
for region in ECN2 ECN1 BRL ECP1 ECP2 FullDet
do
    for build in BH STD CE
    do
	testbase=${base}_${region}_${build}
	rm -rf validation_${testbase}
	rm -rf log_${testbase}_NVU8int_NTH24_cmsswval.txt 
    done
done

rm -rf ${dir}
