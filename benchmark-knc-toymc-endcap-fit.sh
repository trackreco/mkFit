#! /bin/bash

make -j 12

micdir=/nfsmic/${USER}/tmp

for nth in 1 2 4 8 15 30 60 90 120 150 180 210 240
do
    echo "knc toymc" nth=${nth} "FIT (endcap)"
    ssh mic0 ./mkFit-mic  --endcap-test --read --file-name ${micdir}/simtracks_endcap_1kx10k.bin --fit-std-only --num-thr ${nth} >& log_knc_endcap_1kx10k_FIT_NVU16int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8 16
do
    sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
    make clean
    make -j 12

    echo "knc toymc" nvu=${nvu} "FIT (endcap)"
    ssh mic0 ./mkFit-mic --endcap-test --read --file-name ${micdir}/simtracks_endcap_1kx10k.bin --fit-std-only --num-thr 1 >& log_knc_endcap_1kx10k_FIT_NVU${nvu}_NTH1.txt

    sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

make clean
