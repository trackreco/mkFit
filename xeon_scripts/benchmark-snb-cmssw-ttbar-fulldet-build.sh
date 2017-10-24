#! /bin/bash

make -j 12

dir=/data/nfsmic/slava77/samples/2017/pass-4874f28/initialStep/PU70HS/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU/a/
file=memoryFile.fv3.clean.writeAll.recT.072617.bin

base=SNB_CMSSW_TTbar
ttbar=PU70
nevents=240
maxth=24
maxvu=8
exe="./mkFit/mkFit --input-file ${dir}/${file} --cmssw-n2seeds"

for nth in 1 2 4 6 8 12 16 20 24
do
    for bV in "BH bh" "STD std" "CE ce"
    do echo $bV | while read -r bN bO
	do
	    oBase=${base}_${ttbar}_${bN}
	    nproc=$(( ${nevents} / (${maxth}/${nth}) ))
	    bExe="${exe} --build-${bO} --num-thr ${nth} --num-events ${nproc}"

	    echo "${oBase}: Benchmark [nTH:${nth}, nVU:8int]"
	    ${bExe} --dump-for-plots >& log_${oBase}_NVU8int_NTH${nth}.txt

	    if [ "${bN}" == "CE" ]
	    then
		for nev in 1 2 4 8 12
		do
		    if (( ${nev} <= ${nth} ))
		    then
			echo "${oBase}: Benchmark [nTH:${nth}, nVU:8int, nEV:${nev}]"
			${bExe} --silent --num-thr-ev ${nev} >& log_${oBase}_NVU8int_NTH${nth}_NEV${nev}.txt
		    fi
		done
	    fi
	done
    done
done

for nvu in 1 2 4 8
do
    make clean
    make -j 12 USE_INTRINSICS:=-DMPT_SIZE=${nvu}
    for bV in "BH bh" "STD std" "CE ce"
    do echo $bV | while read -r bN bO
	do
	    oBase=${base}_${ttbar}_${bN}
	    nproc=$(( ${nevents} / (${maxvu}/${nvu}) ))
	    bExe="${exe} --build-${bO} --num-thr 1 --num-events ${nproc}"

	    echo "${oBase}: Benchmark [nTH:1, nVU:${nvu}]"
	    ${bExe} --dump-for-plots >& log_${oBase}_NVU${nvu}_NTH1.txt
	done
    done
done

make distclean
