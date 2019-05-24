#! /bin/bash

## input
suite=${1:-"forPR"}
lnxuser=${2:-${USER}}
useLNX=${3:-0}

## In case this is run separately from the main script
source xeon_scripts/common-variables.sh ${suite} ${lnxuser}
source xeon_scripts/init-env.sh

Base_Test="NVU1_NTH1"

if [[ ${useLNX} -eq 0 ]]
then
arch_array=("SNB ${Base_Test}" "SNB NVU8int_NTH24" "KNL ${Base_Test}" "KNL NVU16int_NTH256" "SKL-SP ${Base_Test}" "SKL-SP NVU16int_NTH64")
fi
if [[ ${useLNX} -eq 1 ]]
then
arch_array=("LNX-G ${Base_Test}" "LNX-G NVU16int_NTH64" "LNX-S ${Base_Test}" "LNX-S NVU16int_NTH64")
fi
if [[ ${useLNX} -eq 2 ]]
then
arch_array=("SNB ${Base_Test}" "SNB NVU8int_NTH24" "KNL ${Base_Test}" "KNL NVU16int_NTH256" "SKL-SP ${Base_Test}" "SKL-SP NVU16int_NTH64" "LNX-G ${Base_Test}" "LNX-G NVU16int_NTH64" "LNX-S ${Base_Test}" "LNX-S NVU16int_NTH64")
fi

##### Make plots of track properties (kinematics, nHits, etc) from text files, comparing different machine configurations #####
for build in "${text_builds[@]}"
do echo ${!build} | while read -r bN bO
    do
	echo "Making plots from text files for" ${sample} ":" ${bN}
	
	for archV in "${arch_array[@]}" #"SNB ${Base_Test}" "SNB NVU8int_NTH24" "KNL ${Base_Test}" "KNL NVU16int_NTH256" "SKL-SP ${Base_Test}" "SKL-SP NVU16int_NTH64"
	do echo ${archV} | while read -r archN archO
	    do
		echo "Extracting plots from dump for" ${archN} ${archO}
		python plotting/makePlotsFromDump.py ${archN} ${sample} ${bN} ${archO}
	    done
	done
	
	echo "Making comparison plots from dump for" ${sample} ":" ${bN}
	root -b -q -l plotting/makePlotsFromDump.C\(\"${sample}\",\"${bN}\",\"${suite}\",${useLNX}\)
    done
done
