#! /bin/bash

make -j 12

dir=/data/nfsmic/${USER}/tmp

for nth in 1 2 4 6 8 12 16 20 24
do
    echo "snb toymc" nth=${nth} "FIT (endcap)"
    ./mkFit/mkFit --endcap-test --read --file-name ${dir}/simtracks_endcap_1kx10k.bin --fit-std-only --num-thr ${nth} >& log_snb_endcap_1kx10k_FIT_NVU8int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8
do
    sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
    make clean
    make -j 12

    echo "snb toymc" nvu=${nvu} "FIT (endcap)" 
    ./mkFit/mkFit --endcap-test --read --file-name ${dir}/simtracks_endcap_1kx10k.bin --fit-std-only --num-thr 1 >& log_snb_endcap_1kx10k_FIT_NVU${nvu}_NTH1.txt

    sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

make clean
