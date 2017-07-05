#! /bin/bash

sed -i 's/#WITH_ROOT := yes/WITH_ROOT := yes/g' Makefile.config

make -j 12

dir=/data/nfsmic/scratch/10mu-new

## ECN ##
echo "SNB CMSSW BH (ECN): validation [nTH:24, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${dir}/mu_ecn-1000-10.bin-5 --build-bh  --num-thr 24 >& log_SNB_CMSSW_ECN_BH_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_CMSSW_ECN_BH.root
echo "SNB CMSSW STD (ECN): validation [nTH:24, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${dir}/mu_ecn-1000-10.bin-5 --build-std --num-thr 24 >& log_SNB_CMSSW_ECN_STD_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_CMSSW_ECN_STD.root
echo "SNB CMSSW CE (ECN): validation [nTH:24, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${dir}/mu_ecn-1000-10.bin-5 --build-ce  --num-thr 24 >& log_SNB_CMSSW_ECN_CE_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_CMSSW_ECN_CE.root

## BRL ##
echo "SNB CMSSW BH (BRL): validation [nTH:24, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${dir}/mu_brl-1000-10.bin-5 --build-bh  --num-thr 24 >& log_SNB_CMSSW_BRL_BH_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_CMSSW_BRL_BH.root
echo "SNB CMSSW STD (BRL): validation [nTH:24, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${dir}/mu_brl-1000-10.bin-5 --build-std --num-thr 24 >& log_SNB_CMSSW_BRL_STD_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_CMSSW_BRL_STD.root
echo "SNB CMSSW CE (BRL): validation [nTH:24, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${dir}/mu_brl-1000-10.bin-5 --build-ce  --num-thr 24 >& log_SNB_CMSSW_BRL_CE_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_CMSSW_BRL_CE.root

## ECP ##
echo "SNB CMSSW BH (ECP): validation [nTH:24, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${dir}/mu_ecp-1000-10.bin-5 --build-bh  --num-thr 24 >& log_SNB_CMSSW_ECP_BH_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_CMSSW_ECP_BH.root
echo "SNB CMSSW STD (ECP): validation [nTH:24, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${dir}/mu_ecp-1000-10.bin-5 --build-std --num-thr 24 >& log_SNB_CMSSW_ECP_STD_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_CMSSW_ECP_STD.root
echo "SNB CMSSW CE (ECP): validation [nTH:24, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${dir}/mu_ecp-1000-10.bin-5 --build-ce  --num-thr 24 >& log_SNB_CMSSW_ECP_CE_NVU8int_NTH24_val.txt
mv valtree.root valtree_SNB_CMSSW_ECP_CE.root

sed -i 's/WITH_ROOT := yes/#WITH_ROOT := yes/g' Makefile.config

make clean

for section in ECN BRL ECP
do
    for build in BH STD CE
    do
	root -b -q -l runValidation.C\(\"_SNB_CMSSW_${section}_${build}\"\)
    done
    root -b -q -l makeValidation.C\(\"SNB_CMSSW_${section}\"\)
done

for build in BH STD CE
do
    hadd valtree_SNB_CMSSW_FullDet_${build}.root validation_SNB_CMSSW_ECN_${build}/valtree_SNB_CMSSW_ECN_${build}.root validation_SNB_CMSSW_BRL_${build}/valtree_SNB_CMSSW_BRL_${build}.root validation_SNB_CMSSW_ECP_${build}/valtree_SNB_CMSSW_ECP_${build}.root

    root -b -q -l runValidation.C\(\"_SNB_CMSSW_FullDet_${build}\"\)
done

root -b -q -l makeValidation.C\(\"SNB_CMSSW_FullDet\"\)

make distclean