#! /bin/bash

sed -i 's/\/\/\#define PRINTOUTS_FOR_PLOTS/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make -j 12

micdir=/nfsmic/${USER}/tmp

for nth in 1 2 4 8 15 30 60 90 120 150 180 210 240
do
    echo "knc toymc" nth=${nth} "BH (barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${micdir}/simtracks_barrel_20x10k.bin --build-bh   --num-thr ${nth} >& log_knc_20x10k_BH_NVU16int_NTH${nth}.txt
    echo "knc tomyc" nth=${nth} "COMB (barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${micdir}/simtracks_barrel_20x10k.bin --build-comb --seeds-per-task 32 --num-thr ${nth} --cloner-single-thread >& log_knc_20x10k_COMB_NVU16int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8 16
do
    sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
    make clean
    make -j 12

    echo "knc toymc" nvu=${nvu} "BH (barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${micdir}/simtracks_barrel_20x10k.bin --build-bh  --num-thr 1 >& log_knc_20x10k_BH_NVU${nvu}_NTH1.txt
    echo "knc toymc" nvu=${nvu} "COMB (barrel)"
    ssh mic0 ./mkFit-mic --read --file-name ${micdir}/simtracks_barrel_20x10k.bin --build-comb --seeds-per-task 32 --num-thr 1 --cloner-single-thread >& log_knc_20x10k_COMB_NVU${nvu}_NTH1.txt

    sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

sed -i 's/\#define PRINTOUTS_FOR_PLOTS/\/\/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make clean
