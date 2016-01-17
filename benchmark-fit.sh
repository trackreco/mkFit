#! /bin/bash

sed -i 's/constexpr int nTracks = 20000/constexpr int nTracks = 1000000/g' Config.h 

make clean
make -j 8

for nth in 1 3 7 21
do
echo nth=${nth}
./mkFit/mkFit --fit-std-only --num-thr ${nth} >& log_host_10x1M_FIT_NVU8_NTH${nth}.txt
done

sed -i 's/define MPT_SIZE 8/define MPT_SIZE XX/g' Config.h
for nvu in 1 2 4
do
echo nvu=${nvu}
sed -i "s/define MPT_SIZE XX/define MPT_SIZE ${nvu} \/\/tmp/g" Config.h
make clean
make -j 8
./mkFit/mkFit --fit-std-only --num-thr 1 >& log_host_10x1M_FIT_NVU${nvu}_NTH1.txt
sed -i "s/define MPT_SIZE ${nvu} \/\/tmp/define MPT_SIZE XX/g" Config.h
done

sed -i "s/define MPT_SIZE XX/define MPT_SIZE 8/g" Config.h

sed -i 's/constexpr int nTracks = 1000000/constexpr int nTracks = 20000/g' Config.h 

make clean
make -j 8
