#! /bin/bash

make -j 12 WITH_ROOT=yes

dir=/data/nfsmic/slava77/samples/2017/pass-4874f28/initialStep
file=memoryFile.fv3.recT.072617.bin

ECN2=${dir}/10muEta-24to-17Pt1to10/${file}
ECN1=${dir}/10muEta-175to-055Pt1to10/${file}
BRL=${dir}/10muEtaLT06Pt1to10/${file}
ECP1=${dir}/10muEta055to175Pt1to10/${file}
ECP2=${dir}/10muEta17to24Pt1to10/${file}

base=SNB_CMSSW_10mu

for section in ECN2 ECN1 BRL ECP1 ECP2 
do
    for bV in "BH bh" "STD std" "CE ce" "FV fv"
    do echo $bV | while read -r bN bO
	do
	    oBase=${base}_${section}_${bN}
	    echo "${oBase}: validation [nTH:24, nVU:8]"
	    ./mkFit/mkFit --cmssw-n2seeds --cmssw-val-trkparam --input-file ${!section} --build-${bO} --num-thr 24 >& log_${oBase}_NVU8int_NTH24_cmsswval.txt
	    mv valtree.root valtree_${oBase}.root
	done
    done
done

make clean

for section in ECN2 ECN1 BRL ECP1 ECP2
do
    oBase=${base}_${section}
    for build in BH STD CE FV
    do
	root -b -q -l plotting/runValidation.C\(\"_${oBase}_${build}\",0,1\)
    done
    root -b -q -l plotting/makeValidation.C\(\"${oBase}\",\"\",1\)
done

for build in BH STD CE FV
do
    oBase=${base}
    fBase=valtree_${oBase}
    dBase=validation_${oBase}
    hadd ${fBase}_FullDet_${build}.root `for section in ECN2 ECN1 BRL ECP1 ECP2; do echo -n ${dBase}_${section}_${build}/${fBase}_${section}_${build}.root" "; done`
    root -b -q -l plotting/runValidation.C\(\"_${oBase}_FullDet_${build}\",0,1\)
done
root -b -q -l plotting/makeValidation.C\(\"${oBase}_FullDet\",\"\",1\)

make distclean
