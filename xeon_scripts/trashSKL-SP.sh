#! /bin/bash

useLNX=${1:-0}

# Final cleanup script for benchmarks
if [[ ${useLNX} -eq 3 ]] || [[ ${useLNX} -eq 4 ]] 
then
rm -rf benchmark_knl_dump.txt benchmark_snb_dump.txt
fi
if [[ ${useLNX} -eq 1 ]] || [[ ${useLNX} -eq 2 ]] || [[ ${useLNX} -eq 4 ]]
then
rm -rf benchmark_lnx-g_dump.txt benchmark_lnx-s_dump.txt
fi

rm -rf log_*.txt
rm -rf *.root
rm -rf *.png
rm -rf validation_*
