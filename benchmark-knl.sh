#! /bin/bash

sed -i 's/\/\/\#define PRINTOUTS_FOR_PLOTS/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make -j 64

# dir=/data/nfsmic/${USER}/tmp
# micdir=/nfsmic/${USER}/tmp

# mkdir -p ${dir}
./mkFit/mkFit --write --file-name simtracks_10x20k.bin
# mv simtracks_10x20k.bin ${dir}/

BEXE="./mkFit/mkFit --read --file-name simtracks_10x20k.bin"
# BEXE="numactl -m 1 ./mkFit/mkFit --read --file-name simtracks_10x20k.bin"

for nth in 1 3 7 21 42 63 84 105 126 147 168 189 210 231 252
do
echo "knl" nth=${nth} "BH"
${BEXE} --build-bh  --num-thr ${nth} >& log_knl_10x20k_BH_NVU16int_NTH${nth}.txt
echo "knl" nth=${nth} "STD"
${BEXE} --build-std --num-thr ${nth} >& log_knl_10x20k_ST_NVU16int_NTH${nth}.txt
echo "knl" nth=${nth} "CE"
${BEXE} --build-ce  --num-thr ${nth} >& log_knl_10x20k_CE_NVU16int_NTH${nth}.txt
echo "knl" nth=${nth} "CEST"
${BEXE} --build-ce  --num-thr ${nth} --cloner-single-thread >& log_knl_10x20k_CEST_NVU16int_NTH${nth}.txt
echo "knl" nth=${nth} "TBB"
${BEXE} --build-tbb --seeds-per-task 32 --num-thr ${nth} --cloner-single-thread >& log_knl_10x20k_TBBST_NVU16int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8 16
do
sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
make clean
make -j 64
echo "knl" nvu=${nvu} "BH"
${BEXE} --build-bh  --num-thr 1 >& log_knl_10x20k_BH_NVU${nvu}_NTH1.txt
echo "knl" nvu=${nvu} "STD"
${BEXE} --build-std --num-thr 1 >& log_knl_10x20k_ST_NVU${nvu}_NTH1.txt
echo "knl" nvu=${nvu} "CE"
${BEXE} --build-ce  --num-thr 1 >& log_knl_10x20k_CE_NVU${nvu}_NTH1.txt
echo "knl" nvu=${nvu} "CEST"
${BEXE} --build-ce  --num-thr 1 --cloner-single-thread >& log_knl_10x20k_CEST_NVU${nvu}_NTH1.txt
echo "knl" nvu=${nvu} "TBB"
${BEXE} --build-tbb --seeds-per-task 32 --num-thr 1 --cloner-single-thread >& log_knl_10x20k_TBBST_NVU${nvu}_NTH1.txt
sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

sed -i 's/\#define PRINTOUTS_FOR_PLOTS/\/\/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make clean
make -j 64
