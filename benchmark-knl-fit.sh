#! /bin/bash

sed -i 's/int nTracks = 20000/int nTracks = 1000000/g' Config.cc 

make -j 8

# dir=/data/nfsmic/${USER}/tmp
# micdir=/nfsmic/${USER}/tmp

# mkdir -p ${dir}
#./mkFit/mkFit --write --file-name simtracks_1kx10k.bin
# mv simtracks_10x1M.bin ${dir}/

. data-dir-location.sh

BEXE="./mkFit/mkFit --read --file-name ${dir}/simtracks_barrel_1kxx10k.bin"
# BEXE="numactl -m 1 ./mkFit/mkFit --read --file-name simtracks_10x1M.bin"

LOG_BASE=log_KNL_ToyMC_Barrel

for nth in 1 2 4 8 15 30 60 90 120 150 180 210 240
do
echo "KNL" nth=${nth} "FIT"
${BEXE} --fit-std-only --num-thr ${nth} >& ${LOG_BASE}_FIT_NVU16int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8 16
do
sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
make clean
make -j 8
echo "KNL" nvu=${nvu} "FIT"
${BEXE} --fit-std-only --num-thr 1 >& ${LOG_BASE}_FIT_NVU${nvu}_NTH1.txt
sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

sed -i 's/int nTracks = 1000000/int nTracks = 20000/g' Config.cc 

make clean
make -j 8
