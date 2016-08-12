#! /bin/bash

dir=/nfsmic/cerati/

patch < cmssw-endcap-config.patch

make -j 8

for nth in 1 3 7 21 42 63 84 105 126 147 168 189 210
do
echo "mic cmssw" nth=${nth} "BH (endcap)"
ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-bh  --num-thr ${nth} >& log_mic_endcap_100xTTbarPU35_BH_NVU16int_NTH${nth}.txt
echo "mic cmssw" nth=${nth} "STD (endcap)"
ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-std --num-thr ${nth} >& log_mic_endcap_100xTTbarPU35_ST_NVU16int_NTH${nth}.txt
echo "mic cmssw" nth=${nth} "CE (endcap)"
ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-ce  --num-thr ${nth} >& log_mic_endcap_100xTTbarPU35_CE_NVU16int_NTH${nth}.txt
echo "mic cmssw" nth=${nth} "CEST (endcap)"
ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-ce  --num-thr ${nth} --cloner-single-thread >& log_mic_endcap_100xTTbarPU35_CEST_NVU16int_NTH${nth}.txt
echo "mic cmssw" nth=${nth} "TBB (endcap)"
ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-tbb --seeds-per-task 32 --num-thr ${nth} --cloner-single-thread >& log_mic_endcap_100xTTbarPU35_TBBST_NVU16int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8 16
do
sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
make clean
make -j 8
echo "mic cmssw" nvu=${nvu} "BH (endcap)"
ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-bh  --num-thr 1 >& log_mic_endcap_100xTTbarPU35_BH_NVU${nvu}_NTH1.txt
echo "mic cmssw" nvu=${nvu} "STD (endcap)"
ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-std --num-thr 1 >& log_mic_endcap_100xTTbarPU35_ST_NVU${nvu}_NTH1.txt
echo "mic cmssw" nvu=${nvu} "CE (endcap)"
ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-ce  --num-thr 1 >& log_mic_endcap_100xTTbarPU35_CE_NVU${nvu}_NTH1.txt
echo "mic cmssw" nvu=${nvu} "CEST (endcap)"
ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-ce  --num-thr 1 --cloner-single-thread >& log_mic_endcap_100xTTbarPU35_CEST_NVU${nvu}_NTH1.txt
echo "mic cmssw" nvu=${nvu} "TBB (endcap)"
ssh mic0 ./mkFit-mic --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-tbb --seeds-per-task 32 --num-thr 1 --cloner-single-thread >& log_mic_endcap_100xTTbarPU35_TBBST_NVU${nvu}_NTH1.txt
sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

patch -R < cmssw-endcap-config.patch

make clean
make -j 8
