#! /bin/bash

## In case this is run separately from the main script
[ -z "$ROOTSYS" ] && source ~matevz/root/bin/thisroot.sh
source xeon_scripts/common_variables.sh

Base_Test="NVU1_NTH1"

for build in BH STD CE # FV
do 
    echo "Making plots from text files for" ${physics_sample} ":" ${build}
    
    for archV in "SNB ${Base_Test}" "SNB NVU8int_NTH24" "KNL ${Base_Test}" "KNL NVU16int_NTH256"
    do echo ${archV} | while read -r archN archO
	do
	    echo "Extracting plots from dump for" ${archN} ${archO}
	    python plotting/makePlotsFromDump.py ${archN} ${physics_sample} ${build} ${archO}
	done
    done

    echo "Making comparison plots from dump for" ${physics_sample} ":" ${build}
    root -b -q -l plotting/makePlotsFromDump.C\(\"${physics_sample}\",\"${build}\"\)
done
