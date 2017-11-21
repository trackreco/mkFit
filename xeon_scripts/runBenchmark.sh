#! /bin/bash

# Initialize Benchmarks
[ -z "$ROOTSYS" ] && source ~matevz/root/bin/thisroot.sh
source xeon_scripts/common_variables.sh

##### Launch Tests
echo "Run benchmarking on KNL concurrently with SNB and KNC benchmarks" 
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build-knl.sh >& benchmark_knl_dump.txt &

echo "Run benchmarking on SNB"
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build.sh SNB

echo "Run benchmarking on KNC"
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build.sh KNC

##### Benchmark Plots #####
for archV in "SNB snb" "KNC knc" "KNL knl"
   do echo ${archV} | while read -r archN archO
	do
	echo "Extract benchmarking results for" ${archN}
	python plotting/makeBenchmarkPlots.py ${archN} ${sample}

	echo "Make final plot comparing different build options for" ${archN}
	root -b -q -l plotting/makeBenchmarkPlots.C\(\"${archN}\",\"${sample}\"\)	

	echo "Extract multiple events in flight benchmark results for" ${archN}
	python plotting/makeMEIFBenchmarkPlots.py ${archN} ${sample}

	echo "Make final plot comparing multiple events in flight for" ${archN}
	root -b -q -l plotting/makeMEIFBenchmarkPlots.C\(\"${archN}\",\"${sample}\"\)	
    done
done

##### Plots from Text Files #####
for build in BH STD CE
do 
    echo "Making plots from text files for" ${sample} ":" ${build}
    
    for archV in "SNB NVU8int_NTH24" "KNC NVU16int_NTH240" "KNL NVU16int_NTH256"
    do echo ${archV} | while read -r archN archO
	do
	    echo "Extracting plots from dump for" ${archN} NVU1_NTH1 
	    python plotting/makePlotsFromDump.py ${archN} ${sample} ${build} NVU1_NTH1

	    echo "Extracting plots from dump for" ${archN} ${archO}
	    python plotting/makePlotsFromDump.py ${archN} ${sample} ${build} ${archO}
	done
    done

    echo "Making comparison plots from dump for" ${sample} ":" ${build}
    root -b -q -l plotting/makePlotsFromDump.C\(\"${sample}\",\"${build}\"\)
done

##### Validation tests #####
echo "Running ROOT based validation"
./val_scripts/validation-snb-cmssw-ttbar-fulldet-build-allval.sh

##### Final cleanup #####
make distclean
