#! /bin/bash

[ -z "$ROOTSYS" ] && source ~matevz/root/bin/thisroot.sh

sample="CMSSW_TTbar_PU70"

##### Benchmark Tests #####
for archV in "SNB snb" "KNC knc" "KNL knl"
    do echo $archV | while read -r archN archO
	do
	echo "Run benchmarking on" ${archN}
	./xeon_scripts/benchmark-${archO}-cmssw-ttbar-fulldet-build.sh
	
	echo "Extract benchmarking results for" ${archN}
	python plotting/makeBenchmarkPlots.py ${archN} ${sample}

	echo "Make final plot comparing different build options for" ${archN}
	root -b -q -l plotting/makeBenchmarkPlots.C\(\"${archN}\",\"${sample}\"\)	

	echo "Extract multiple events in flight benchmark results for" ${archN}
	python plotting/makeMEIFBenchmarkPlots.py ${archN} ${sample}

	echo "Make final plot comparing mulitple events in flight for" ${archN}
	root -b -q -l plotting/makeMEIFBenchmarkPlots.C\(\"${archN}\",\"${sample}\"\)	
    done
done

##### nHits plots #####
for build in BH STD CE
do 
    echo "Making nHits plots for" ${sample} ":" ${build}
    
    for archV in "SNB NVU8int_NTH24" "KNC NVUY16int_NTH240" "KNL NVU16int_NTH256"
    do echo $archV | while read -r archN archO
	do
	    echo "Extracting nHits for" ${archN} NVU1_NTH1 
	    python plotting/makePlotsFromDump.py ${archN} ${sample} ${build} NVU1_NTH1

	    echo "Extracting nHits for" ${archN} ${archO}
	    python plotting/makePlotsFromDump.py ${archN} ${sample} ${build} ${archO}
	done
    done

    echo "Making final plot comparing nHits for" ${sample} ":" ${build}
    root -b -q -l plotting/makePlotsFromDump.C\(\"${sample}\",\"${build}\"\)
done

##### Validation tests #####
echo "Running ROOT based validation"
./val_scripts/validation-snb-cmssw-ttbar-fulldet-build-allval.sh

make distclean
