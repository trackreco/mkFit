#! /bin/bash

sed -i 's/#WITH_ROOT := yes/WITH_ROOT := yes/g' Makefile.config

make -j 12

dir=/data/nfsmic/kmcdermo/toymc

echo "snb toymc BH (barrel): validation [nTH:24, nVU:8]"
./mkFit/mkFit --normal-val --read --file-name ${dir}/simtracks_barrel_10x10k_val.bin --build-bh  --num-thr 24 >& log_snb_10x10k_BH_NVU8int_NTH24_val.txt
mv valtree.root valtree_BH.root
echo "snb toymc STD (barrel): validation [nTH:24, nVU:8]"
./mkFit/mkFit --normal-val --read --file-name ${dir}/simtracks_barrel_10x10k_val.bin --build-std --num-thr 24 >& log_snb_10x10k_STD_NVU8int_NTH24_val.txt
mv valtree.root valtree_STD.root
echo "snb toymc CE (barrel): validation [nTH:24, nVU:8]"
./mkFit/mkFit --normal-val --read --file-name ${dir}/simtracks_barrel_10x10k_val.bin --build-ce  --num-thr 24 --cloner-single-thread >& log_snb_10x10k_CE_NVU8int_NTH24_val.txt
mv valtree.root valtree_CE.root

sed -i 's/WITH_ROOT := yes/#WITH_ROOT := yes/g' Makefile.config

make clean

