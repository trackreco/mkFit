#!/bin/bash

#################################################################################
##                                    README!                                  ##
##                                                                             ##
## Stress test script to run on phi3, testing different thread/MEIF combos     ##
## with different instruction sets using clone engine building and n2-seeding. ##
##                                                                             ##
## 120 events per physical core are processed (need even number for ST tests). ##
## Can vary nevents, thread/MEIF combos, input file, seeds, building algo.     ##
##                                                                             ##
## Command line inputs are which platform to stress (ben_arch) no_turbo ON/OFF ##
## (no_turbo), the min time per test (duration), and the time between each     ##
## test (sleep_time).                                                          ##
##                                                                             ##    
## Output file lists stress test time per event processed per physical core.   ##
#################################################################################

########################
## Source Environment ##
########################

source /opt/rh/devtoolset-7/enable
source /opt/intel/bin/compilervars.sh intel64
source stress-test-config.sh

###################
## Configuration ##
###################

## Command line inputs
ben_arch=${1} # SNB (phi1), KNL (phi2), SKL-SP (phi3)
no_turbo=${2:-1} # turbo on/off (default is turbo OFF)
duration=${3:-1800} # min time spent for each test [s]
sleep_time=${4:-300} # sleep time between tests [s]

## platform specific settings
if [[ "${ben_arch}" == "SNB" ]]
then
    mOpt="-j 12"
    dir=/data2/nfsmic/slava77/samples
    maxcore=12
    declare -a instruction_sets=(SSE3 AVX)
    declare -a thread_combo_arr=("1 1" "6 6" "12 6" "12 12" "24 6" "24 12" "24 24")
    declare -a njob_arr=("12" "24")
elif [[ "${ben_arch}" == "KNL" ]]
then
    mOpt="-j 64"
    dir=/data1/work/slava77/samples
    maxcore=64
    declare -a instruction_sets=(SSE3 AVX AX2 AVX512)
    declare -a thread_combo_arr=("1 1" "32 32" "64 32" "64 64" "128 32" "128 64" "128 128" "256 32" "256 64" "256 128" "256 256")
    declare -a njob_arr=("32" "64" "128" "256")
elif [[ "${ben_arch}" == "SKL-SP" ]]
then
    mOpt="-j 32"
    dir=/data2/slava77/samples
    maxcore=32
    declare -a instruction_sets=(SSE3 AVX AX2 AVX512)
    declare -a thread_combo_arr=("1 1" "16 16" "32 16" "32 32" "48 16" "48 32" "64 16" "64 32" "64 64")
    declare -a njob_arr=("32" "64")
else 
    echo "${ben_arch} is not a valid architecture! Exiting..."
    exit
fi

## Common file setup
subdir=2017/pass-c93773a/initialStep/PU70HS/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU
file=memoryFile.fv3.clean.writeAll.CCC1620.recT.082418-25daeda.bin
nevents=120

## Common mkFit options
seeds="--cmssw-n2seeds"
algo="--ce"
opts="--silent"
exe="./mkFit/mkFit --input-file ${dir}/${subdir}/${file} ${seeds} ${algo} ${opts}"

## Output options
base_outname="mkfit_stress_tests_${ben_arch}"
output_file="${base_outname}_results.${ext}"

##################
## Set no_turbo ##
##################

echo ${no_turbo} | sudo /usr/bin/tee /sys/devices/system/cpu/intel_pstate/no_turbo >/dev/null 2>&1

###############
## Run tests ##
###############

## loop instruction sets (i.e. build minimally)
for instruction_set in "${instruction_sets[@]}"
do
    ## compile once, using settings for the given instruction set
    make distclean
    make "${mOpt}" "${!instruction_set}"
    
    ## run thread combo tests (nThreads, nEventsInFlight)
    for thread_combo in "${thread_combo_arr[@]}"
    do echo "${thread_combo}" | while read -r nth nev
	do
	    ## compute total number of events to process
	    ncore=$( GetNCore() "${nth}" "${maxcore}" ) 
	    nproc=$(( ${nevents} * ${ncore} ))

	    ## print out which test is being performed
	    test_label="${instruction_set}_${nth_label}${nth}_${nev_label}${nev}"
	    echo "Running stress test for: ${test_label}..."

	    ## full executable
	    base_exe="${exe} --num-thr ${nth} --num-thr-ev ${nev}"

	    ## output file
	    tmp_output_file="${base_outname}_${test_label}.${tmp_ext}"
	    
	    ## execute test and pipe time to output file: https://stackoverflow.com/a/2409214
	    { time MkFitLoop "${duration}" "${base_exe}" "${nproc}" "1" > /dev/null 2>&1 ; } 2> "${tmp_output_file}"

	    ## pause to let machine cool down between each test
	    sleep "${sleep_time}"

	    ## add other info about test to tmp file
	    AppendTmpFile "${tmp_output_file}" "${ncore}" "${nproc}" "${nloop}"
	done # end loop over reading thread combo
    done # end loop over thread combos

    ## run special test of N jobs, single thread each
    for njob in "${njob_arr[@]}"
    do
	## compute total number of events to process
	ncore=$( GetNCore() "${njob}" "${maxcore}" ) 
	nproc=$(( ${nevents} * ${ncore} ))

	## print out which test is being performed
	test_label="${instruction_set}_${njob_label}${njob}"
	echo "Running stress test for: ${test_label}..."

	## base executable
	base_exe="${exe} --num-thr 1 --num-thr-ev 1"

	## output file
	tmp_output_file="${base_outname}_${test_label}.${tmp_ext}"
	    
	## execute test and pipe time to output file: https://stackoverflow.com/a/2409214
	{ time MkFitLoop "${duration}" "${base_exe}" "${nproc}" "${njob}" > /dev/null 2>&1 ; } 2> "${tmp_output_file}"

        ## add other info about test to tmp file
	AppendTmpFile "${tmp_output_file}" "${ncore}" "${nproc}" "${nloop}"

	## pause to let machine cool down between each test
	sleep "${sleep_time}"
    done # end loop over njob for single thread

done # end loop over instruction set

#######################
## Make Final Output ##
#######################

## final output file
> "${output_file}"

## loop over all output files, and append results to single file
for thread_combo in "${thread_combo_arr[@]}"
do echo "${thread_combo}" | while read -r nth nev
    do
	## get test label, print it
	test_label="${instruction_set}_${nth_label}${nth}_${nev_label}${nev}"
	echo "Computing time for: ${test_label}"
	
        ## get tmp output file name
	tmp_output_file="${base_outname}_${test_label}.${tmp_ext}"

	## dump into output file
	DumpIntoFile "${tmp_output_file}" "${output_file}"
    done
done

## loop over single thread njob tests, and append to single file
for njob in "${njob_arr[@]}"
do
    ## get test label, print it
    test_label="${instruction_set}_${njob_label}${njob}"
    echo "Computing time for: ${test_label}"

    ## get tmp output file name
    tmp_output_file="${base_outname}_${test_label}.${tmp_ext}"

    ## dump into output file
    DumpIntoFile "${tmp_output_file}" "${output_file}"
done

####################
## Reset no_turbo ##
####################

## default state of machine is turbo OFF
echo 1 | sudo /usr/bin/tee /sys/devices/system/cpu/intel_pstate/no_turbo >/dev/null 2>&1 

###################
## Final Message ##
###################

echo "Finished stress test!"
