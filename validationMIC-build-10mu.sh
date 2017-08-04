#! /bin/bash


[ -e "$BIN_DATA_PATH" ] || BIN_DATA_PATH=/store/disk00/slava77/analysis/CMSSW_9_1_0_pre1-tkNtuple/run1000/2017/pass-4874f28/initialStep
ECN2=${BIN_DATA_PATH}/10muEta-24to-17Pt1to10/memoryFile.fv3.recT.072617.bin
ECN1=${BIN_DATA_PATH}/10muEta-175to-055Pt1to10/memoryFile.fv3.recT.072617.bin
BRL=${BIN_DATA_PATH}/10muEtaLT06Pt1to10/memoryFile.fv3.recT.072617.bin
ECP1=${BIN_DATA_PATH}/10muEta055to175Pt1to10/memoryFile.fv3.recT.072617.bin
ECP2=${BIN_DATA_PATH}/10muEta17to24Pt1to10/memoryFile.fv3.recT.072617.bin

runValidation(){
    for sV in "sim " "see --cmssw-seeds"; do echo $sV | while read -r sN sO; do
	    for section in ECN2 ECN1 BRL ECP1 ECP2; do
	        for bV in "BH bh" "STD std" "CE ce"; do echo $bV | while read -r bN bO; do
		        oBase=${base}_${sN}_${section}_${bN}
		        nTH=8
		        echo "${oBase}: validation [nTH:${nTH}, nVU:8]"
		        ./mkFit/mkFit --geom CMS-2017 --root-val --read --file-name ${!section} --build-${bO} ${sO} --num-thr ${nTH} >& log_${oBase}_NVU8int_NTH${nTH}_val.txt
		        mv valtree.root valtree_${oBase}.root
		    done
	        done
	    done
        done
    done
    
    make clean
    
    for opt in sim see
    do
        for section in ECN2 ECN1 BRL ECP1 ECP2
        do
	    oBase=${base}_${opt}_${section}
	    for build in BH STD CE
	    do
	        root -b -q -l runValidation.C+\(\"_${oBase}_${build}\"\)
	    done
	    root -b -q -l makeValidation.C+\(\"${oBase}\"\)
        done
        
        for build in BH STD CE
        do
	    oBase=${base}_${opt}
	    fBase=valtree_${oBase}
	    dBase=validation_${oBase}
	    hadd ${fBase}_FullDet_${build}.root `for section in ECN2 ECN1 BRL ECP1 ECP2; do echo -n ${dBase}_${section}_${build}/${fBase}_${section}_${build}.root" "; done`
            
	    root -b -q -l runValidation.C+\(\"_${oBase}_FullDet_${build}\"\)
        done
        
        root -b -q -l makeValidation.C+\(\"${oBase}_FullDet\"\)
        
    done
    
    make distclean
}

#cleanup first
make clean
make distclean

export base=SNB_CMSSW_10mu
echo Run default with base = ${base}
make -j 12 WITH_ROOT=yes
#this does make clean inside
runValidation

export base=SNB_CMSSW_10mu_cleanSeed
echo Run CLEAN_SEEDS with base = ${base}
make -j 12 WITH_ROOT=yes CPPUSERFLAGS+=-DCLEAN_SEEDS CXXUSERFLAGS+=-DCLEAN_SEEDS
#this does make clean inside
runValidation

unset base
