#! /bin/bash

sed -i 's/constexpr int nTracks = 20000/constexpr int nTracks = 1000000/g' Config.h 

make clean
make -j 8

for nth in 1 3 7 21
do
echo nth=${nth}
./mkFit/mkFit --fit-std-only --num-thr ${nth} >& log_host_10x1M_FIT_NVU8int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8
do
sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
make clean
make -j 8
echo nvu=${nvu}
./mkFit/mkFit --fit-std-only --num-thr 1 >& log_host_10x1M_FIT_NVU${nvu}_NTH1.txt
sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

sed -i 's/constexpr int nTracks = 1000000/constexpr int nTracks = 20000/g' Config.h 

make clean
make -j 8
