#! /bin/bash


dir=/data/nfsmic/cerati/

patch < cmssw-config.patch

make -j 8

for nth in 1 3 7 21
do
echo "host cmssw" nth=${nth} "BH"
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-bh  --num-thr ${nth} >& log_host_100xTTbarPU35_BH_NVU8int_NTH${nth}.txt
echo "host cmssw" nth=${nth} "STD"
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-std --num-thr ${nth} >& log_host_100xTTbarPU35_ST_NVU8int_NTH${nth}.txt
echo "host cmssw" nth=${nth} "CE"
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-ce  --num-thr ${nth} >& log_host_100xTTbarPU35_CE_NVU8int_NTH${nth}.txt
echo "host cmssw" nth=${nth} "CEST"
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-ce  --num-thr ${nth} --cloner-single-thread >& log_host_100xTTbarPU35_CEST_NVU8int_NTH${nth}.txt
echo "host cmssw" nth=${nth} "TBB"
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-tbb --num-thr ${nth} --cloner-single-thread >& log_host_100xTTbarPU35_TBBST_NVU8int_NTH${nth}.txt
done
for nth in 10 12 14 16
do
echo "host cmssw" nth=${nth} "TBB"
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-tbb --num-thr ${nth} --cloner-single-thread >& log_host_100xTTbarPU35_TBBST_NVU8int_NTH${nth}.txt
echo "host cmssw" nth=${nth} "BH"
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-bh  --num-thr ${nth} >& log_host_100xTTbarPU35_BH_NVU8int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8
do
sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
make clean
make -j 8
echo "host cmssw" nvu=${nvu} "BH"
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-bh  --num-thr 1 >& log_host_100xTTbarPU35_BH_NVU${nvu}_NTH1.txt
echo "host cmssw" nvu=${nvu} "STD"
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-std --num-thr 1 >& log_host_100xTTbarPU35_ST_NVU${nvu}_NTH1.txt
echo "host cmssw" nvu=${nvu} "CE"
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-ce  --num-thr 1 >& log_host_100xTTbarPU35_CE_NVU${nvu}_NTH1.txt
echo "host cmssw" nvu=${nvu} "CEST"
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-ce  --num-thr 1 --cloner-single-thread >& log_host_100xTTbarPU35_CEST_NVU${nvu}_NTH1.txt
echo "host cmssw" nvu=${nvu} "TBB"
./mkFit/mkFit --read --file-name ${dir}/cmssw_100xTTbarPU35_polar_split_mock_noFWD.bin --cms-geom --cmssw-seeds --build-tbb --num-thr 1 --cloner-single-thread >& log_host_100xTTbarPU35_TBBST_NVU${nvu}_NTH1.txt
sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

patch -R < cmssw-config.patch

make clean
make -j 8
