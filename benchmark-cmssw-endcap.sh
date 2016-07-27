#! /bin/bash


dir=/data/nfsmic/cerati/

patch < cmssw-endcap-config.patch

make -j 8

for nth in 1 3 7 21
do
echo nth=${nth}
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-bh  --num-thr ${nth} >& log_host_endcap_100xTTbarPU35_BH_NVU8int_NTH${nth}.txt
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-std --num-thr ${nth} >& log_host_endcap_100xTTbarPU35_ST_NVU8int_NTH${nth}.txt
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-ce  --num-thr ${nth} >& log_host_endcap_100xTTbarPU35_CE_NVU8int_NTH${nth}.txt
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-ce  --num-thr ${nth} --cloner-single-thread >& log_host_endcap_100xTTbarPU35_CEST_NVU8int_NTH${nth}.txt
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-tbb --num-thr ${nth} --cloner-single-thread >& log_host_endcap_100xTTbarPU35_TBBST_NVU8int_NTH${nth}.txt
done
for nth in 10 12 14 16
do
echo nth=${nth}
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-tbb --num-thr ${nth} --cloner-single-thread >& log_host_endcap_100xTTbarPU35_TBBST_NVU8int_NTH${nth}.txt
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-bh  --num-thr ${nth} >& log_host_endcap_100xTTbarPU35_BH_NVU8int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8
do
sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
make clean
make -j 8
echo nvu=${nvu}
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-bh  --num-thr 1 >& log_host_endcap_100xTTbarPU35_BH_NVU${nvu}_NTH1.txt
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-std --num-thr 1 >& log_host_endcap_100xTTbarPU35_ST_NVU${nvu}_NTH1.txt
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-ce  --num-thr 1 >& log_host_endcap_100xTTbarPU35_CE_NVU${nvu}_NTH1.txt
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-ce  --num-thr 1 --cloner-single-thread >& log_host_endcap_100xTTbarPU35_CEST_NVU${nvu}_NTH1.txt
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_endcapEta1p7.bin --endcap-test --cms-geom --cmssw-seeds --build-tbb --num-thr 1 --cloner-single-thread >& log_host_endcap_100xTTbarPU35_TBBST_NVU${nvu}_NTH1.txt
sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

patch -R < cmssw-endcap-config.patch

make clean
make -j 8
