#! /bin/bash

sed -i 's/\/\/\#define PRINTOUTS_FOR_PLOTS/\#define PRINTOUTS_FOR_PLOTS/g' mkFit/MkBuilder.cc

make clean
make -j 8

for nth in 1 3 7 21 #42 63 84 105 126 147 168 189 210
do
echo nth=${nth}
ssh mic0 ./mkFit-mic --build-bh --num-thr ${nth} >& log_mic_10x20k_BH_NVU16_NTH${nth}.txt
ssh mic0 ./mkFit-mic --build-ce --num-thr ${nth} >& log_mic_10x20k_CE_NVU16_NTH${nth}.txt
ssh mic0 ./mkFit-mic --build-ce --num-thr ${nth} --cloner-single-thread >& log_mic_10x20k_CEST_NVU16_NTH${nth}.txt
done

#sed -i 's/define MPT_SIZE 16/define MPT_SIZE XX/g' Config.h
#for nvu in 1 2 4 8
#do
#sed -i "s/define MPT_SIZE XX/define MPT_SIZE ${nvu} \/\/tmp/g" Config.h
#make clean
#make -j 8
#echo nvu=${nvu}
#ssh mic0 ./mkFit-mic --build-bh --num-thr 1 >& log_mic_10x20k_BH_NVU${nvu}_NTH1.txt
#ssh mic0 ./mkFit-mic --build-ce --num-thr 1 >& log_mic_10x20k_CE_NVU${nvu}_NTH1.txt
#ssh mic0 ./mkFit-mic --build-ce --num-thr 1 --cloner-single-thread >& log_mic_10x20k_CEST_NVU${nvu}_NTH1.txt
#sed -i "s/define MPT_SIZE ${nvu} \/\/tmp/define MPT_SIZE XX/g" Config.h
#done

#sed -i "s/define MPT_SIZE XX/define MPT_SIZE 16/g" Config.h

sed -i 's/\#define PRINTOUTS_FOR_PLOTS/\/\/\#define PRINTOUTS_FOR_PLOTS/g' mkFit/MkBuilder.cc

make clean
make -j 8
