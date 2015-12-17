#! /bin/bash

sed -i 's/\/\/\#define PRINTOUTS_FOR_PLOTS/\#define PRINTOUTS_FOR_PLOTS/g' mkFit/MkBuilder.cc

make clean
make -j 8

for nth in 1 3 7
do
./mkFit/mkFit --build-bh --num-thr ${nth} >& log_10x20k_BH_NVU8_NTH${nth}.txt
./mkFit/mkFit --build-ce --num-thr ${nth} >& log_10x20k_CE_NVU8_NTH${nth}.txt
./mkFit/mkFit --build-ce --num-thr ${nth} --cloner-single-thread >& log_10x20k_CEST_NVU8_NTH${nth}.txt
done

sed -i 's/define MPT_SIZE 8/define MPT_SIZE XX/g' Config.h
for nvu in 2 4 8
do
sed -i "s/define MPT_SIZE XX/define MPT_SIZE ${nvu}/g" Config.h
make clean
make -j 8
./mkFit/mkFit --build-bh --num-thr 1 >& log_10x20k_BH_NVU${nvu}_NTH1.txt
./mkFit/mkFit --build-ce --num-thr 1 >& log_10x20k_CE_NVU${nvu}_NTH1.txt
./mkFit/mkFit --build-ce --num-thr 1 --cloner-single-thread >& log_10x20k_CEST_NVU${nvu}_NTH1.txt
sed -i "s/define MPT_SIZE ${nvu}/define MPT_SIZE XX/g" Config.h
done

sed -i "s/define MPT_SIZE XX/define MPT_SIZE 8/g" Config.h

sed -i 's/\#define PRINTOUTS_FOR_PLOTS/\/\/\#define PRINTOUTS_FOR_PLOTS/g' mkFit/MkBuilder.cc

make clean
make -j 8
