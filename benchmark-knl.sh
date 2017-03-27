#! /bin/bash

sed -i 's/\/\/\#define PRINTOUTS_FOR_PLOTS/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make -j 64

# dir=/data/nfsmic/${USER}/tmp
# micdir=/nfsmic/${USER}/tmp

# mkdir -p ${dir}
#./mkFit/mkFit --write --file-name simtracks_20x10k.bin
# mv simtracks_10x20k.bin ${dir}/

BEXE="./mkFit/mkFit --read --file-name simtracks_20x10k.bin"
# BEXE="numactl -m 1 ./mkFit/mkFit --read --file-name simtracks_10x20k.bin"

LOG_BASE=log_KNL_ToyMC_Barrel

for nth in 1 2 4 8 15 30 60 90 120 150 180 210 240
do
echo "KNL" nth=${nth} "BH"
${BEXE} --build-bh  --num-thr ${nth} >& ${LOG_BASE}_BH_NVU16int_NTH${nth}.txt
echo "KNL" nth=${nth} "STD"
${BEXE} --build-std --num-thr ${nth} >& ${LOG_BASE}_ST_NVU16int_NTH${nth}.txt
echo "KNL" nth=${nth} "CE"
${BEXE} --build-ce  --num-thr ${nth} >& ${LOG_BASE}_CE_NVU16int_NTH${nth}.txt
echo "KNL" nth=${nth} "CEST"
${BEXE} --build-ce  --num-thr ${nth} --cloner-single-thread >& ${LOG_BASE}_CEST_NVU16int_NTH${nth}.txt
echo "KNL" nth=${nth} "TBB"
${BEXE} --build-tbb --seeds-per-task 32 --num-thr ${nth} --cloner-single-thread >& ${LOG_BASE}_TBBST_NVU16int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8 16
do
sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
make clean
make -j 64
echo "KNL" nvu=${nvu} "BH"
${BEXE} --build-bh  --num-thr 1 >& ${LOG_BASE}_BH_NVU${nvu}_NTH1.txt
echo "KNL" nvu=${nvu} "STD"
${BEXE} --build-std --num-thr 1 >& ${LOG_BASE}_ST_NVU${nvu}_NTH1.txt
echo "KNL" nvu=${nvu} "CE"
${BEXE} --build-ce  --num-thr 1 >& ${LOG_BASE}_CE_NVU${nvu}_NTH1.txt
echo "KNL" nvu=${nvu} "CEST"
${BEXE} --build-ce  --num-thr 1 --cloner-single-thread >& ${LOG_BASE}_CEST_NVU${nvu}_NTH1.txt
echo "KNL" nvu=${nvu} "TBB"
${BEXE} --build-tbb --seeds-per-task 32 --num-thr 1 --cloner-single-thread >& ${LOG_BASE}_TBBST_NVU${nvu}_NTH1.txt
sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

sed -i 's/\#define PRINTOUTS_FOR_PLOTS/\/\/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make clean
make -j 64
