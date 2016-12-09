#! /bin/bash

sed -i 's/\/\/\#define PRINTOUTS_FOR_PLOTS/\#define PRINTOUTS_FOR_PLOTS/g' Config.h
patch < cmssw-config.patch

make -j 12

dir=/nfsmic/cerati/

for nth in 1 2 4 8 15 30 60 90 120 150 180 210 240
do
    echo "KNC CMSSW" nth=${nth} "BH (Barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-bh  --num-thr ${nth} >& log_KNC_CMSSW_Barrel_BH_NVU16int_NTH${nth}.txt
    echo "KNC CMSSW" nth=${nth} "STD (Barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-std --seeds-per-task 32 --num-thr ${nth} >& log_KNC_CMSSW_Barrel_STD_NVU16int_NTH${nth}.txt
    echo "KNC CMSSW" nth=${nth} "CE (Barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-ce  --seeds-per-task 32 --num-thr ${nth} >& log_KNC_CMSSW_Barrel_CE_NVU16int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8 16
do
    sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
    make clean
    make -j 12
   
    echo "KNC CMSSW" nvu=${nvu} "BH (Barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-bh  --num-thr 1 >& log_KNC_CMSSW_Barrel_BH_NVU${nvu}_NTH1.txt
    echo "KNC CMSSW" nvu=${nvu} "STD (Barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-std --seeds-per-task 32 --num-thr 1 >& log_KNC_CMSSW_Barrel_STD_NVU${nvu}_NTH1.txt
    echo "KNC CMSSW" nvu=${nvu} "CE (Barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-ce  --seeds-per-task 32 --num-thr 1 >& log_KNC_CMSSW_Barrel_CE_NVU${nvu}_NTH1.txt

    sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

patch -R < cmssw-config.patch
sed -i 's/\#define PRINTOUTS_FOR_PLOTS/\/\/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make clean

