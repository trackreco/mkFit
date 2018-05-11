#! /bin/bash

make -j 32 WITH_ROOT=yes

dir=/data2/slava77/samples/2017/pass-4874f28/initialStep
file=memoryFile.fv3.clean.writeAll.recT.072617.bin

NoPU=10024.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017_GenSimFullINPUT+DigiFull_2017+RecoFull_2017+ALCAFull_2017+HARVESTFull_2017
PU35=10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU
PU70=PU70/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU

base=SKL-SP_CMSSW_TTbar

for ttbar in NoPU PU35 PU70 
do
    for sV in "SimSeed --cmssw-simseeds" "CMSSeed --cmssw-n2seeds"
    do echo $sV | while read -r sN sO
	do
	    for bV in "BH bh" "STD std" "CE ce" "FV fv"
	    do echo $bV | while read -r bN bO
		do
		    oBase=${base}_${ttbar}_${sN}_${bN}
		    echo "${oBase}: validation [nTH:32, nVU:32]"
		    ./mkFit/mkFit ${sO} --sim-val --input-file ${dir}/${!ttbar}/${file} --build-${bO} --num-thr 32 >& log_${oBase}_NVU32int_NTH32_val.txt
		    mv valtree.root valtree_${oBase}.root
		done
	    done
	done
    done
done

make clean

for ttbar in NoPU PU35 PU70 
do
    for seed in SimSeed CMSSeed
    do
	oBase=${base}_${ttbar}_${seed}
	for build in BH STD CE FV
	do
	    root -b -q -l plotting/runValidation.C\(\"_${oBase}_${build}\"\)
	done
	root -b -q -l plotting/makeValidation.C\(\"${oBase}\"\)
    done
done

make distclean
