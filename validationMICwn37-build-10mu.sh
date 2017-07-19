#! /bin/bash

sed -i 's/#WITH_ROOT := yes/WITH_ROOT := yes/g' Makefile.config

make -j 12

ecn2=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep/10muEta-24to-17Pt1to10/memoryFile.fv3.recT.071817.bin
ecn1=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep/10muEta-175to-055Pt1to10/memoryFile.fv3.recT.071817.bin
brl=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep/10muEtaLT06Pt1to10/memoryFile.fv3.recT.071817.bin
ecp1=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep/10muEta055to175Pt1to10/memoryFile.fv3.recT.071817.bin
ecp2=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep/10muEta17to24Pt1to10/memoryFile.fv3.recT.071817.bin


## ECN2 ##
echo "SNB CMSSW BH (ECN2): validation [nTH:1, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${ecn2} --build-bh  --num-thr 8 >& log_SNB_CMSSW_ECN2_BH_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_ECN2_BH.root
echo "SNB CMSSW STD (ECN2): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${ecn2} --build-std --num-thr 8 >& log_SNB_CMSSW_ECN2_STD_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_ECN2_STD.root
echo "SNB CMSSW CE (ECN2): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${ecn2} --build-ce  --num-thr 8 >& log_SNB_CMSSW_ECN2_CE_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_ECN2_CE.root

## ECN1 ##
echo "SNB CMSSW BH (ECN1): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${ecn1} --build-bh  --num-thr 8 >& log_SNB_CMSSW_ECN1_BH_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_ECN1_BH.root
echo "SNB CMSSW STD (ECN1): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${ecn1} --build-std --num-thr 8 >& log_SNB_CMSSW_ECN1_STD_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_ECN1_STD.root
echo "SNB CMSSW CE (ECN1): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${ecn1} --build-ce  --num-thr 8 >& log_SNB_CMSSW_ECN1_CE_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_ECN1_CE.root

## BRL ##
echo "SNB CMSSW BH (BRL): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${brl} --build-bh  --num-thr 8 >& log_SNB_CMSSW_BRL_BH_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_BRL_BH.root
echo "SNB CMSSW STD (BRL): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${brl} --build-std --num-thr 8 >& log_SNB_CMSSW_BRL_STD_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_BRL_STD.root
echo "SNB CMSSW CE (BRL): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${brl} --build-ce  --num-thr 8 >& log_SNB_CMSSW_BRL_CE_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_BRL_CE.root

## ECP1 ##
echo "SNB CMSSW BH (ECP1): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${ecp1} --build-bh  --num-thr 8 >& log_SNB_CMSSW_ECP1_BH_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_ECP1_BH.root
echo "SNB CMSSW STD (ECP1): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${ecp1} --build-std --num-thr 8 >& log_SNB_CMSSW_ECP1_STD_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_ECP1_STD.root
echo "SNB CMSSW CE (ECP1): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${ecp1} --build-ce  --num-thr 8 >& log_SNB_CMSSW_ECP1_CE_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_ECP1_CE.root

## ECP2 ##
echo "SNB CMSSW BH (ECP2): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${ecp2} --build-bh  --num-thr 8 >& log_SNB_CMSSW_ECP2_BH_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_ECP2_BH.root
echo "SNB CMSSW STD (ECP2): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${ecp2} --build-std --num-thr 8 >& log_SNB_CMSSW_ECP2_STD_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_ECP2_STD.root
echo "SNB CMSSW CE (ECP2): validation [nTH:8, nVU:8]"
./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${ecp2} --build-ce  --num-thr 8 >& log_SNB_CMSSW_ECP2_CE_NVU8int_NTH8_val.txt
mv valtree.root valtree_SNB_CMSSW_ECP2_CE.root



sed -i 's/WITH_ROOT := yes/#WITH_ROOT := yes/g' Makefile.config

make clean

for section in ECN2 ECN1 BRL ECP1 ECP2
do
    for build in BH STD CE
    do
	root -b -q -l runValidation.C+\(\"_SNB_CMSSW_${section}_${build}\"\)
    done
    root -b -q -l makeValidation.C+\(\"SNB_CMSSW_${section}\"\)
done

for build in BH STD CE
do
    hadd valtree_SNB_CMSSW_FullDet_${build}.root validation_SNB_CMSSW_ECN2_${build}/valtree_SNB_CMSSW_ECN2_${build}.root validation_SNB_CMSSW_ECN1_${build}/valtree_SNB_CMSSW_ECN1_${build}.root validation_SNB_CMSSW_BRL_${build}/valtree_SNB_CMSSW_BRL_${build}.root validation_SNB_CMSSW_ECP1_${build}/valtree_SNB_CMSSW_ECP1_${build}.root validation_SNB_CMSSW_ECP2_${build}/valtree_SNB_CMSSW_ECP2_${build}.root

    root -b -q -l runValidation.C+\(\"_SNB_CMSSW_FullDet_${build}\"\)
done

root -b -q -l makeValidation.C+\(\"SNB_CMSSW_FullDet\"\)

make distclean