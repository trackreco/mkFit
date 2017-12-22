#! /bin/bash

make -j 12 WITH_ROOT=yes

dir=/data/nfsmic/scratch/toymc
file=simtracks_fulldet_400x2p5k_val.bin

base=SNB_ToyMC_FullDet

for bV in "BH bh" "STD std" "CE ce" "FV fv"
do echo $bV | while read -r bN bO
    do
	oBase=${base}_${bN}
	echo "${oBase}: validation [nTH:24, nVU:8]"
	./mkFit/mkFit --root-val --read-simtrack-states --seed-input sim --input-file ${dir}/${file} --build-${bO} --num-thr 24 >& log_${oBase}_NVU8int_NTH24_val.txt
	mv valtree.root valtree_${oBase}.root
    done
done

for build in BH STD CE FV
do
    root -b -q -l plotting/runValidation.C\(\"_SNB_ToyMC_FullDet_${build}\",1\)
done
root -b -q -l plotting/makeValidation.C\(\"SNB_ToyMC_FullDet\"\)

make clean
