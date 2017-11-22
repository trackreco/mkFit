#! /bin/bash

[ -e "$BIN_DATA_PATH" ] || BIN_DATA_PATH=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep
fin=${BIN_DATA_PATH}/PU70/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU/memoryFile.fv3.clean.writeAll.recT.072617.bin

runValidation()
{
    for sV in "sim --cmssw-simseeds" "see --cmssw-stdseeds"; do echo $sV | while read -r sN sO; do
	    if [ "${1}" == "1" ]; then
		sO="--cmssw-n2seeds"
	    fi
            for bV in "BH bh" "STD std" "CE ce" "FV fv"; do echo $bV | while read -r bN bO; do
		    oBase=${base}_${sN}_${bN}
		    nTH=8
		    echo "${oBase}: validation [nTH:${nTH}, nVU:8]"
		    ./mkFit/mkFit --root-val --input-file ${fin} --build-${bO} ${sO} --num-thr ${nTH} >& log_${oBase}_NVU8int_NTH${nTH}_val.txt
		    mv valtree.root valtree_${oBase}.root
                done
            done
        done
    done
        
    for opt in sim see
    do
        oBase=${base}_${opt}
        for build in BH STD CE FV
        do
	    root -b -q -l plotting/runValidation.C+\(\"_${oBase}_${build}\"\)
        done
        root -b -q -l plotting/makeValidation.C+\(\"${oBase}\"\)
    done
}

#cleanup first
make clean
make distclean
make -j 12 WITH_ROOT=yes

export base=SNB_CMSSW_PU70_clean
echo Run default build with base = ${base}
runValidation 0

export base=SNB_CMSSW_PU70_clean_cleanSeed
echo Run CLEAN_SEEDS with base = ${base}
runValidation 1

make distclean

unset base
