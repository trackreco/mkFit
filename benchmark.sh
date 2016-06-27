#! /bin/bash

sed -i 's/\/\/\#define PRINTOUTS_FOR_PLOTS/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make -j 8

dir=/data/nfsmic/${USER}/tmp

mkdir -p ${dir}
./mkFit/mkFit --write --file-name simtracks_10x20k.bin
mv simtracks_10x20k.bin ${dir}/

for nth in 1 3 7 21
do
echo "host" nth=${nth} "BH"
./mkFit/mkFit --read --file-name ${dir}/simtracks_10x20k.bin --build-bh  --num-thr ${nth} >& log_host_10x20k_BH_NVU8int_NTH${nth}.txt
echo "host" nth=${nth} "STD"
./mkFit/mkFit --read --file-name ${dir}/simtracks_10x20k.bin --build-std --num-thr ${nth} >& log_host_10x20k_ST_NVU8int_NTH${nth}.txt
echo "host" nth=${nth} "CE"
./mkFit/mkFit --read --file-name ${dir}/simtracks_10x20k.bin --build-ce  --num-thr ${nth} >& log_host_10x20k_CE_NVU8int_NTH${nth}.txt
echo "host" nth=${nth} "CEST"
./mkFit/mkFit --read --file-name ${dir}/simtracks_10x20k.bin --build-ce  --num-thr ${nth} --cloner-single-thread >& log_host_10x20k_CEST_NVU8int_NTH${nth}.txt
echo "host" nth=${nth} "TBB"
./mkFit/mkFit --read --file-name ${dir}/simtracks_10x20k.bin --build-tbb --num-thr ${nth} --cloner-single-thread >& log_host_10x20k_TBBST_NVU8int_NTH${nth}.txt
done
for nth in 10 12 14 16
do
echo "host" nth=${nth} "TBB"
./mkFit/mkFit --read --file-name ${dir}/simtracks_10x20k.bin --build-tbb --num-thr ${nth} --cloner-single-thread >& log_host_10x20k_TBBST_NVU8int_NTH${nth}.txt
echo "host" nth=${nth} "BH"
./mkFit/mkFit --read --file-name ${dir}/simtracks_10x20k.bin --build-bh  --num-thr ${nth} >& log_host_10x20k_BH_NVU8int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8
do
sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
make clean
make -j 8
echo "host" nvu=${nvu} "BH"
./mkFit/mkFit --read --file-name ${dir}/simtracks_10x20k.bin --build-bh  --num-thr 1 >& log_host_10x20k_BH_NVU${nvu}_NTH1.txt
echo "host" nvu=${nvu} "STD"
./mkFit/mkFit --read --file-name ${dir}/simtracks_10x20k.bin --build-std --num-thr 1 >& log_host_10x20k_ST_NVU${nvu}_NTH1.txt
echo "host" nvu=${nvu} "CE"
./mkFit/mkFit --read --file-name ${dir}/simtracks_10x20k.bin --build-ce  --num-thr 1 >& log_host_10x20k_CE_NVU${nvu}_NTH1.txt
echo "host" nvu=${nvu} "CEST"
./mkFit/mkFit --read --file-name ${dir}/simtracks_10x20k.bin --build-ce  --num-thr 1 --cloner-single-thread >& log_host_10x20k_CEST_NVU${nvu}_NTH1.txt
echo "host" nvu=${nvu} "TBB"
./mkFit/mkFit --read --file-name ${dir}/simtracks_10x20k.bin --build-tbb --num-thr 1 --cloner-single-thread >& log_host_10x20k_TBBST_NVU${nvu}_NTH1.txt
sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

sed -i 's/\#define PRINTOUTS_FOR_PLOTS/\/\/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make clean
make -j 8
