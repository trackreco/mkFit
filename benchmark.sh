#! /bin/bash

sed -i 's/\/\/\#define PRINTOUTS_FOR_PLOTS/\#define PRINTOUTS_FOR_PLOTS/g' mkFit/MkBuilder.cc

make clean
make -j 8

for nth in 1 3 7 21
do
echo nth=${nth}
./mkFit/mkFit --build-bh  --num-thr ${nth} >& log_host_10x20k_BH_NVU8_NTH${nth}.txt
./mkFit/mkFit --build-std --num-thr ${nth} >& log_host_10x20k_ST_NVU8_NTH${nth}.txt
./mkFit/mkFit --build-ce  --num-thr ${nth} >& log_host_10x20k_CE_NVU8_NTH${nth}.txt
./mkFit/mkFit --build-ce  --num-thr ${nth} --cloner-single-thread >& log_host_10x20k_CEST_NVU8_NTH${nth}.txt
done

sed -i 's/define MPT_SIZE 8/define MPT_SIZE XX/g' Config.h
for nvu in 1 2 4
do
echo nvu=${nvu}
sed -i "s/define MPT_SIZE XX/define MPT_SIZE ${nvu} \/\/tmp/g" Config.h
make clean
make -j 8
./mkFit/mkFit --build-bh  --num-thr 1 >& log_host_10x20k_BH_NVU${nvu}_NTH1.txt
./mkFit/mkFit --build-std --num-thr 1 >& log_host_10x20k_ST_NVU${nvu}_NTH1.txt
./mkFit/mkFit --build-ce  --num-thr 1 >& log_host_10x20k_CE_NVU${nvu}_NTH1.txt
./mkFit/mkFit --build-ce  --num-thr 1 --cloner-single-thread >& log_host_10x20k_CEST_NVU${nvu}_NTH1.txt
sed -i "s/define MPT_SIZE ${nvu} \/\/tmp/define MPT_SIZE XX/g" Config.h
done

sed -i "s/define MPT_SIZE XX/define MPT_SIZE 8/g" Config.h

sed -i 's/\#define PRINTOUTS_FOR_PLOTS/\/\/\#define PRINTOUTS_FOR_PLOTS/g' mkFit/MkBuilder.cc

make clean
make -j 8
