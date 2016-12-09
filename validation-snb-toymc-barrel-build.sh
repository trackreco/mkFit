#! /bin/bash

sed -i 's/#WITH_ROOT := yes/WITH_ROOT := yes/g' Makefile.config

make -j 12

dir=/data/nfsmic/kmcdermo/toymc

echo "SNB ToyMC BH (Barrel): validation [nTH:24, nVU:8]"
./mkFit/mkFit --normal-val --read --file-name ${dir}/simtracks_barrel_10x10k_val.bin --build-bh  --num-thr 24 >& log_SNB_ToyMC_Barrel_BH_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_ToyMC_Barrel_BH.root
echo "SNB ToyMC STD (Barrel): validation [nTH:24, nVU:8]"
./mkFit/mkFit --normal-val --read --file-name ${dir}/simtracks_barrel_10x10k_val.bin --build-std --num-thr 24 >& log_SNB_ToyMC_Barrel_STD_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_ToyMC_Barrel_STD.root
echo "SNB ToyMC CE (Barrel): validation [nTH:24, nVU:8]"
./mkFit/mkFit --normal-val --read --file-name ${dir}/simtracks_barrel_10x10k_val.bin --build-ce  --num-thr 24 >& log_SNB_ToyMC_Barrel_CE_NVU8int_NTH24_val.txt
mv valtree.root valtree_snb_ToyMC_Barrel_CE.root

sed -i 's/WITH_ROOT := yes/#WITH_ROOT := yes/g' Makefile.config

make clean

