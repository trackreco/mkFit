#! /bin/bash

make -j 12

dir=/data/nfsmic/slava77/samples/2017/pass-4874f28/initialStep/PU70HS/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU/a/
file=memoryFile.fv3.clean.writeAll.recT.072617.bin

base=SNB_CMSSW_TTbar
ttbar=PU70
nevents=240
minth=1
maxth=24
minvu=1
maxvu=8
nevdump=50
dump=DumpForPlots
exe="./mkFit/mkFit --input-file ${dir}/${file} --cmssw-n2seeds"

for nth in 1 2 4 6 8 12 16 20 24
do
    for bV in "BH bh" "STD std" "CE ce"
    do echo $bV | while read -r bN bO
	do
	    oBase=${base}_${ttbar}_${bN}
	    nproc=$(( ${nevents} / (${maxth}/${nth}) ))
	    bExe="${exe} --build-${bO} --num-thr ${nth}"

	    echo "${oBase}: Benchmark [nTH:${nth}, nVU:${maxvu}int]"
	    ${bExe} --num-events ${nproc} >& log_${oBase}_NVU${maxvu}int_NTH${nth}.txt

	    if [ "${bN}" == "CE" ]
	    then
		for nev in 1 2 4 8 12
		do
		    if (( ${nev} <= ${nth} ))
		    then
			echo "${oBase}: Benchmark [nTH:${nth}, nVU:${maxvu}int, nEV:${nev}]"
			${bExe} --silent --num-thr-ev ${nev} --num-events ${nproc} >& log_${oBase}_NVU${maxvu}int_NTH${nth}_NEV${nev}.txt
		    fi
		done
	    fi

	    if (( ${nth} == ${maxth} ))
	    then
		echo "${oBase}: Text dump for plots [nTH:${nth}, nVU:${maxvu}int]"
		${bExe} --dump-for-plots --num-events ${nevdump} >& log_${oBase}_NVU${maxvu}int_NTH${nth}_${dump}.txt
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
	    bExe="${exe} --build-${bO} --num-thr ${minth}"

	    echo "${oBase}: Benchmark [nTH:${minth}, nVU:${nvu}]"
	    ${bExe} --num-events ${nproc} >& log_${oBase}_NVU${nvu}_NTH1.txt

	    if (( ${nvu} == ${minvu} ))
	    then
		echo "${oBase}: Text dump for plots [nTH:${minth}, nVU:${nvu}]"
		${bExe} --dump-for-plots --num-events ${nevdump} >& log_${oBase}_NVU${nvu}_NTH${minth}_${dump}.txt
	    fi
	done
    done
done

make distclean
