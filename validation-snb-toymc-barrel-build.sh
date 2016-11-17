#! /bin/bash

sed -i 's/#WITH_ROOT := yes/WITH_ROOT := yes/g' Makefile.config

make -j 12

dir=/data/nfsmic/${USER}/tmp

echo "snb toymc BH (barrel): validation [nTH:24, nVU:8]"
./mkFit/mkFit --normal-val --read --file-name ${dir}/simtracks_barrel_10x10k_val.bin --build-bh   --num-thr 24 >& log_snb_10x10k_BH_NVU8int_NTH24_val.txt
mv valtree.root valtree_BH.root
echo "snb toymc COMB (barrel): validation [nTH:24, nVU:8]"
./mkFit/mkFit --normal-val --read --file-name ${dir}/simtracks_barrel_10x10k_val.bin --build-comb --num-thr 24 --cloner-single-thread >& log_snb_10x10k_COMB_NVU8int_NTH24_val.txt
mv valtree.root valtree_COMB.root

sed -i 's/WITH_ROOT := yes/#WITH_ROOT := yes/g' Makefile.config

make clean

