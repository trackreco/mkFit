#! /bin/bash

# Main output dir name
dir=${1:-benchmarks}

# First collect all plots and text files into common dir
echo "Moving plots and text files locally to ${dir}"
./web/collectBenchmarks.sh ${dir}

# Now copy to lxplus
echo "Moving plots and text files remotely to ${whost}"
./web/copyAndSendToLXPLUS.sh

# Final cleanup of directory
echo "Removing local files"
./xeon_scripts/trashSNB.sh
rm -rf ${dir}
