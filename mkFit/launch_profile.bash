#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-socket=1
#SBATCH --gres=gpu:1
#SBATCH -t 0:30:00

##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user= ml15@princeton.edu


module load cudatoolkit

exec="./mkFit"
date1=$(date +"%s")

echo `which nvprof`
echo $exec


#nvprof --export-profile timeline.nvprof ./a.out
#nvprof --analysis-metrics -o anaylsis.nvprof ./a.out 
#nvprof --metrics achieved_occupancy,executed_ipc -o metrics.nvprof ./a.out  
nvprof --metrics all -o metrics512.nvprof $exec
#nvprof --query-events

date2=$(date +"%s")
diff=$(($date2-$date1))

echo "--------------------------------------------"
echo "Time elapsed: "
echo "$(($diff / 60)):$(($diff % 60)) mimutes"
echo "--------------------------------------------"
