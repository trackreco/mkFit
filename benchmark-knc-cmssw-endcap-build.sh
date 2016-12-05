#! /bin/bash

sed -i 's/\/\/\#define PRINTOUTS_FOR_PLOTS/\#define PRINTOUTS_FOR_PLOTS/g' Config.h
patch < cmssw-endcap-config.patch

make -j 12

dir=/nfsmic/cerati/

for nth in 1 2 4 8 15 30 60 90 120 150 180 210 240
do
    echo "knc cmssw" nth=${nth} "BH (endcap)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-bh  --num-thr ${nth} >& log_mic_endcap_100xTTbarPU35_BH_NVU16int_NTH${nth}.txt
    echo "knc cmssw" nth=${nth} "STD (endcap)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-std --seeds-per-task 32 --num-thr ${nth} >& log_mic_endcap_100xTTbarPU35_STD_NVU16int_NTH${nth}.txt
    echo "knc cmssw" nth=${nth} "CE (endcap)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-ce  --seeds-per-task 32 --num-thr ${nth} --cloner-single-thread >& log_mic_endcap_100xTTbarPU35_CE_NVU16int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8 16
do
    sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
    make clean
    make -j 12

    echo "knc cmssw" nvu=${nvu} "BH (endcap)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-bh  --num-thr 1 >& log_mic_endcap_100xTTbarPU35_BH_NVU${nvu}_NTH1.txt
    echo "knc cmssw" nvu=${nvu} "STD (endcap)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-std --seeds-per-task 32 --num-thr 1 >& log_mic_endcap_100xTTbarPU35_STD_NVU${nvu}_NTH1.txt
    echo "knc cmssw" nvu=${nvu} "CE (endcap)"
    ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-ce  --seeds-per-task 32 --num-thr 1 --cloner-single-thread >& log_mic_endcap_100xTTbarPU35_CE_NVU${nvu}_NTH1.txt

    sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

patch -R < cmssw-endcap-config.patch
sed -i 's/\#define PRINTOUTS_FOR_PLOTS/\/\/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make clean
