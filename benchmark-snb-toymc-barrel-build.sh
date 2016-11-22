#! /bin/bash

sed -i 's/\/\/\#define PRINTOUTS_FOR_PLOTS/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make -j 12

dir=/data/nfsmic/kmcdermo/toymc

for nth in 1 2 4 6 8 12 16 20 24
do
    echo "snb toymc" nth=${nth} "BH (barrel)"
    ./mkFit/mkFit --read --file-name ${dir}/simtracks_barrel_20x10k.bin --build-bh   --num-thr ${nth} >& log_snb_20x10k_BH_NVU8int_NTH${nth}.txt
    echo "snb toymc" nth=${nth} "COMB (barrel)"
    ./mkFit/mkFit --read --file-name ${dir}/simtracks_barrel_20x10k.bin --build-comb --num-thr ${nth} --cloner-single-thread >& log_snb_20x10k_COMB_NVU8int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8
do
    sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
    make clean
    make -j 12
    
    echo "snb toymc" nvu=${nvu} "BH (barrel)"
    ./mkFit/mkFit --read --file-name ${dir}/simtracks_barrel_20x10k.bin --build-bh   --num-thr 1 >& log_snb_20x10k_BH_NVU${nvu}_NTH1.txt
    echo "snb toymc" nvu=${nvu} "COMB (barrel)"
    ./mkFit/mkFit --read --file-name ${dir}/simtracks_barrel_20x10k.bin --build-comb --num-thr 1 --cloner-single-thread >& log_snb_20x10k_COMB_NVU${nvu}_NTH1.txt

    sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

sed -i 's/\#define PRINTOUTS_FOR_PLOTS/\/\/\#define PRINTOUTS_FOR_PLOTS/g' Config.h

make clean

