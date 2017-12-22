#! /bin/bash

[ -e "$BIN_DATA_PATH" ] || BIN_DATA_PATH=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep
fin=${BIN_DATA_PATH}/PU70/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU/memoryFile.fv3.clean.writeAll.recT.072617.bin

runBenchmark()
{
#    for sV in "sim --cmssw-simseeds" "see --cmssw-stdseeds"; do echo $sV | while read -r sN sO; do
    for sV in "see --cmssw-stdseeds"; do echo $sV | while read -r sN sO; do
            if [ "${1}" == "1" ]; then
                sO="--cmssw-n2seeds"
            fi
            for bV in "BH bh" "STD std" "CE ce" "FV fv"; do echo $bV | while read -r bN bO; do
		    oBase=${base}_${sN}_${bN}
		    for nTH in 1 4 8 16 32; do
		        echo "${oBase}: benchmark [nTH:${nTH}, nVU:8]"
		        time ./mkFit/mkFit --input-file ${fin} --build-${bO} ${sO} --num-thr ${nTH} >& log_${oBase}_NVU8int_NTH${nTH}_benchmark.txt
		    done
                done
            done
        done
    done
}

#cleanup first
make clean
make distclean

make -j 12
export base=SNB_CMSSW_PU70_clean
echo Run default build with base = ${base}
runBenchmark 0

export base=SNB_CMSSW_PU70_clean_cleanSeed
echo Run CLEAN_SEEDS build with base = ${base}
runBenchmark 1
make clean
make distclean


make -j 12 CPPUSERFLAGS+="-march=native -mtune=native" CXXUSERFLAGS+="-march=native -mtune=native"
export base=SNB_CMSSW_PU70_clean_native
echo Run native build with base = ${base}
runBenchmark 0

export base=SNB_CMSSW_PU70_clean_native_cleanSeed
echo Run CLEAN_SEEDS build with base = ${base}
runBenchmark 1
make clean
make distclean


ECN2=${BIN_DATA_PATH}/10muEta-24to-17Pt1to10/memoryFile.fv3.recT.072617.bin
ECN1=${BIN_DATA_PATH}/10muEta-175to-055Pt1to10/memoryFile.fv3.recT.072617.bin
BRL=${BIN_DATA_PATH}/10muEtaLT06Pt1to10/memoryFile.fv3.recT.072617.bin
ECP1=${BIN_DATA_PATH}/10muEta055to175Pt1to10/memoryFile.fv3.recT.072617.bin
ECP2=${BIN_DATA_PATH}/10muEta17to24Pt1to10/memoryFile.fv3.recT.072617.bin

runBenchmarkSections()
{
    for sV in "sim --cmssw-seeds" "see --cmssw-stdseeds"; do echo $sV | while read -r sN sO; do
            if [ "${1}" == "1" ]; then
                sO="--cmssw-n2seeds"
            fi
            for section in ECN2 ECN1 BRL ECP1 ECP2; do
                for bV in "BH bh" "STD std" "CE ce" "FV fv"; do echo $bV | while read -r bN bO; do
                        oBase=${base}_${sN}_${section}_${bN}
                        nTH=8
                        echo "${oBase}: benchmark [nTH:${nTH}, nVU:8]"
                        time ./mkFit/mkFit --input-file ${!section} --build-${bO} ${sO} --num-thr ${nTH} >& log_${oBase}_NVU8int_NTH${nTH}_benchmark.txt
                    done
                done
            done
        done
    done

}

#this part has a pretty limited value due to the tiny load in the muon samples
make -j 12
export base=SNB_CMSSW_10mu
echo Run default build with base = ${base}
runBenchmarkSections 1

make clean
make distclean


unset base

