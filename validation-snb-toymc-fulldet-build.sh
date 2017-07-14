#! /bin/bash

sed -i 's/#WITH_ROOT := yes/WITH_ROOT := yes/g' Makefile.config

make -j 12

dir=/data/nfsmic/scratch/toymc

echo "SNB ToyMC BH (FullDet): validation [nTH:24, nVU:8]"
./mkFit/mkFit --root-val --read --file-name ${dir}/simtracks_fulldet_400x2p5k_val.bin --build-bh  --num-thr 24 >& log_SNB_ToyMC_FullDet_BH_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_ToyMC_FullDet_BH.root
echo "SNB ToyMC STD (FullDet): validation [nTH:24, nVU:8]"
./mkFit/mkFit --root-val --read --file-name ${dir}/simtracks_fulldet_400x2p5k_val.bin --build-std --num-thr 24 >& log_SNB_ToyMC_FullDet_STD_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_ToyMC_FullDet_STD.root
echo "SNB ToyMC CE (FullDet): validation [nTH:24, nVU:8]"
./mkFit/mkFit --root-val --read --file-name ${dir}/simtracks_fulldet_400x2p5k_val.bin --build-ce  --num-thr 24 >& log_SNB_ToyMC_FullDet_CE_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_ToyMC_FullDet_CE.root

sed -i 's/WITH_ROOT := yes/#WITH_ROOT := yes/g' Makefile.config

for build in BH STD CE
do
    root -b -q -l runValidation.C\(\"_SNB_ToyMC_FullDet_${build}\"\)
done
root -b -q -l makeValidation.C\(\"SNB_ToyMC_FullDet\"\)

make clean
