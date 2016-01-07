#! /bin/bash

sed -i 's/\/\/\#define PRINTOUTS_FOR_PLOTS/\#define PRINTOUTS_FOR_PLOTS/g' mkFit/MkBuilder.cc

make clean
make -j 8

for nth in 1 3 7 21 #42 63 84 105 126 147 168 189 210
do
echo nth=${nth}
ssh mic0 ./mkFit-mic --build-bh --num-thr ${nth} >& log_mic_10x20k_BH_NVU16int_NTH${nth}.txt
ssh mic0 ./mkFit-mic --build-ce --num-thr ${nth} >& log_mic_10x20k_CE_NVU16int_NTH${nth}.txt
ssh mic0 ./mkFit-mic --build-ce --num-thr ${nth} --cloner-single-thread >& log_mic_10x20k_CEST_NVU16int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8 16
do
sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
make clean
make -j 8
echo nvu=${nvu}
ssh mic0 ./mkFit-mic --build-bh --num-thr 1 >& log_mic_10x20k_BH_NVU${nvu}_NTH1.txt
ssh mic0 ./mkFit-mic --build-ce --num-thr 1 >& log_mic_10x20k_CE_NVU${nvu}_NTH1.txt
ssh mic0 ./mkFit-mic --build-ce --num-thr 1 --cloner-single-thread >& log_mic_10x20k_CEST_NVU${nvu}_NTH1.txt
sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

sed -i 's/\#define PRINTOUTS_FOR_PLOTS/\/\/\#define PRINTOUTS_FOR_PLOTS/g' mkFit/MkBuilder.cc

make clean
make -j 8
