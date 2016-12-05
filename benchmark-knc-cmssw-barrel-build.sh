#! /bin/bash

sed -i 's/\/\/\#define PRINTOUTS_FOR_PLOTS/\#define PRINTOUTS_FOR_PLOTS/g' Config.h
patch < cmssw-config.patch

make -j 12

dir=/nfsmic/cerati/

for nth in 1 2 4 8 15 30 60 90 120 150 180 210 240
do
    echo "knc cmssw" nth=${nth} "BH (barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-bh  --num-thr ${nth} >& log_knc_100xTTbarPU35_BH_NVU16int_NTH${nth}.txt
    echo "knc cmssw" nth=${nth} "STD (barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-std --seeds-per-task 32 --num-thr ${nth} >& log_knc_100xTTbarPU35_STD_NVU16int_NTH${nth}.txt
    echo "knc cmssw" nth=${nth} "CE (barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-ce  --seeds-per-task 32 --num-thr ${nth} --cloner-single-thread >& log_knc_100xTTbarPU35_CE_NVU16int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8 16
do
    sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
    make clean
    make -j 12
    
    echo "knc cmssw" nvu=${nvu} "BH (barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-bh  --num-thr 1 >& log_knc_100xTTbarPU35_BH_NVU${nvu}_NTH1.txt
    echo "knc cmssw" nvu=${nvu} "STD (barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-std --seeds-per-task 32 --num-thr 1 >& log_knc_100xTTbarPU35_STD_NVU${nvu}_NTH1.txt
    echo "knc cmssw" nvu=${nvu} "CE (barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-ce  --seeds-per-task 32 --num-thr 1 --cloner-single-thread >& log_knc_100xTTbarPU35_CE_NVU${nvu}_NTH1.txt

    sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

patch -R < cmssw-config.patch
sed -i 's/\#define PRINTOUTS_FOR_PLOTS/\/\/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make clean

