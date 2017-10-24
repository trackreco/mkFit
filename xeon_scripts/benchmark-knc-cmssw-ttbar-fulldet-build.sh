#! /bin/bash

mKNC="KNC_BUILD:=1"
make -j 12 ${mKNC}

dir=/nfsmic/slava77/samples/2017/pass-4874f28/initialStep/PU70HS/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU/a/
file=memoryFile.fv3.clean.writeAll.recT.072617.bin

base=KNC_CMSSW_TTbar
ttbar=PU70
nevents=2400
maxth=240
maxvu=16
exe="ssh mic0 ./mkFit-mic --input-file ${dir}/${file} --cmssw-n2seeds"

for nth in 1 2 4 8 15 30 60 90 120 150 180 210 240
do
    for bV in "BH bh" "STD std" "CE ce"
    do echo $bV | while read -r bN bO
	do
	    oBase=${base}_${ttbar}_${bN}
	    nproc=$(( ${nevents} / (${maxth}/${nth}) ))
	    bExe="${exe} --build-${bO} --num-thr ${nth} --num-events ${nproc}"

	    echo "${oBase}: Benchmark [nTH:${nth}, nVU:16int]"
	    ${bExe} --dump-for-plots >& log_${oBase}_NVU16int_NTH${nth}.txt

	    if [ "${bN}" == "CE" ]
	    then
		for nev in 1 2 4 8 16 32 64 128
		do
		    if (( ${nev} <= ${nth} ))
		    then
			echo "${oBase}: Benchmark [nTH:${nth}, nVU:16int, nEV:${nev}]"
			${bExe} --silent --num-thr-ev ${nev} >& log_${oBase}_NVU16int_NTH${nth}_NEV${nev}.txt
		    fi
		done
	    fi
	done
    done
done

for nvu in 1 2 4 8 16
do
    make clean ${mKNC}
    make -j 12 ${mKNC} USE_INTRINSICS:=-DMPT_SIZE=${nvu}
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

make distclean ${mKNC}
