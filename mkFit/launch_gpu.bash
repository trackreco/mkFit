#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-socket=1
#SBATCH --gres=gpu:1
#SBATCH -t 0:15:00

##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=ml15@princeton.edu

#./sgemm_cublas 
#./sgemm_linear_shared
#./sgemm_mplex_shared_streamed
./mkFit
#./sgemm_mplex_shared
