#! /bin/bash

## input
suite=${1:-"forPR"}

## In case this is run separately from the main script
source xeon_scripts/common-variables.sh ${suite}
source xeon_scripts/init-env.sh

Base_Test="NVU1_NTH1"

##### Make plots of track properties (kinematics, nHits, etc) from text files, comparing different machine configurations #####
for build in "${text_builds[@]}"
do echo ${!build} | while read -r bN bO
    do
	echo "Making plots from text files for" ${sample} ":" ${bN}
	
	for archV in "SNB ${Base_Test}" "SNB NVU8int_NTH24" "KNL ${Base_Test}" "KNL NVU16int_NTH256" "SKL-SP ${Base_Test}" "SKL-SP NVU16int_NTH64"
	do echo ${archV} | while read -r archN archO
	    do
		echo "Extracting plots from dump for" ${archN} ${archO}
		python plotting/makePlotsFromDump.py ${archN} ${sample} ${bN} ${archO}
	    done
	done
	
	echo "Making comparison plots from dump for" ${sample} ":" ${bN}
	root -b -q -l plotting/makePlotsFromDump.C\(\"${sample}\",\"${bN}\",\"${suite}\"\)
    done
done
