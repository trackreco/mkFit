#! /bin/bash

make -j 12 WITH_ROOT=yes CPPUSERFLAGS+=-DCLEAN_SEEDS CXXUSERFLAGS+=-DCLEAN_SEEDS

dir=/data/nfsmic/slava77/samples/2017/pass-4874f28/initialStep/PU70/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU
file=memoryFile.fv3.clean.writeAll.recT.072617.bin

base=SNB_CMSSW_PU70

for sV in "SimSeed " "CMSSeed --cmssw-seeds"
do echo $sV | while read -r sN sO
    do
	for bV in "BH bh" "STD std" "CE ce"
	do echo $bV | while read -r bN bO
	    do
		oBase=${base}_${sN}_${bN}
		echo "${oBase}: validation [nTH:24, nVU:8]"
		./mkFit/mkFit ${sO} --geom CMS-2017 --root-val --read --file-name ${dir}/${file} --build-${bO} --num-thr 24 >& log_${oBase}_NVU8int_NTH24_val.txt
		mv valtree.root valtree_${oBase}.root
	    done
	done
    done
done

make clean

for seed in SimSeed CMSSeed
do
    oBase=${base}_${seed}
    for build in BH STD CE
    do
	root -b -q -l runValidation.C\(\"_${oBase}_${build}\"\)
    done
    root -b -q -l makeValidation.C\(\"${oBase}\"\)
done

make distclean
