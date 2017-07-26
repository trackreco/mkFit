#! /bin/bash

make -j 12 WITH_ROOT=yes

ECN2=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep/10muEta-24to-17Pt1to10/memoryFile.fv3.recT.071817.bin
ECN1=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep/10muEta-175to-055Pt1to10/memoryFile.fv3.recT.071817.bin
BRL=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep/10muEtaLT06Pt1to10/memoryFile.fv3.recT.071817.bin
ECP1=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep/10muEta055to175Pt1to10/memoryFile.fv3.recT.071817.bin
ECP2=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep/10muEta17to24Pt1to10/memoryFile.fv3.recT.071817.bin

for section in ECN2 ECN1 BRL ECP1 ECP2; do
    echo "SNB CMSSW_sim BH (${section}): validation [nTH:1, nVU:8]"
    ./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${!section} --build-bh  --num-thr 8 >& log_SNB_CMSSW_sim_${section}_BH_NVU8int_NTH8_val.txt
    mv valtree.root valtree_SNB_CMSSW_sim_${section}_BH.root
    echo "SNB CMSSW_sim STD (${section}): validation [nTH:8, nVU:8]"
    ./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${!section} --build-std --num-thr 8 >& log_SNB_CMSSW_sim_${section}_STD_NVU8int_NTH8_val.txt
    mv valtree.root valtree_SNB_CMSSW_sim_${section}_STD.root
    echo "SNB CMSSW_sim CE (${section}): validation [nTH:8, nVU:8]"
    ./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${!section} --build-ce  --num-thr 8 >& log_SNB_CMSSW_sim_${section}_CE_NVU8int_NTH8_val.txt
    mv valtree.root valtree_SNB_CMSSW_sim_${section}_CE.root
done
for section in ECN2 ECN1 BRL ECP1 ECP2; do
    echo "SNB CMSSW_see BH (${section}): validation [nTH:1, nVU:8]"
    ./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${!section} --build-bh --cmssw-seeds --num-thr 8 >& log_SNB_CMSSW_see_${section}_BH_NVU8int_NTH8_val.txt
    mv valtree.root valtree_SNB_CMSSW_see_${section}_BH.root
    echo "SNB CMSSW_see STD (${section}): validation [nTH:8, nVU:8]"
    ./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${!section} --build-std --cmssw-seeds --num-thr 8 >& log_SNB_CMSSW_see_${section}_STD_NVU8int_NTH8_val.txt
    mv valtree.root valtree_SNB_CMSSW_see_${section}_STD.root
    echo "SNB CMSSW_see CE (${section}): validation [nTH:8, nVU:8]"
    ./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${!section} --build-ce --cmssw-seeds --num-thr 8 >& log_SNB_CMSSW_see_${section}_CE_NVU8int_NTH8_val.txt
    mv valtree.root valtree_SNB_CMSSW_see_${section}_CE.root
done

make clean

for opt in sim see
do
    for section in ECN2 ECN1 BRL ECP1 ECP2
    do
	for build in BH STD CE
	do
	    root -b -q -l runValidation.C+\(\"_SNB_CMSSW_${opt}_${section}_${build}\"\)
	done
	root -b -q -l makeValidation.C+\(\"SNB_CMSSW_${opt}_${section}\"\)
    done

    for build in BH STD CE
    do
	hadd valtree_SNB_CMSSW_${opt}_FullDet_${build}.root validation_SNB_CMSSW_${opt}_ECN2_${build}/valtree_SNB_CMSSW_${opt}_ECN2_${build}.root validation_SNB_CMSSW_${opt}_ECN1_${build}/valtree_SNB_CMSSW_${opt}_ECN1_${build}.root validation_SNB_CMSSW_${opt}_BRL_${build}/valtree_SNB_CMSSW_${opt}_BRL_${build}.root validation_SNB_CMSSW_${opt}_ECP1_${build}/valtree_SNB_CMSSW_${opt}_ECP1_${build}.root validation_SNB_CMSSW_${opt}_ECP2_${build}/valtree_SNB_CMSSW_${opt}_ECP2_${build}.root

	root -b -q -l runValidation.C+\(\"_SNB_CMSSW_${opt}_FullDet_${build}\"\)
    done

    root -b -q -l makeValidation.C+\(\"SNB_CMSSW_${opt}_FullDet\"\)

done

make distclean

