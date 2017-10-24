#! /bin/bash

# To be sourced where needed

host=`hostname`

if [[ $host == phi2.t2.* ]]; then
  dir=/data/scratch/toymc
  n_sim_thr=128
elif [[ $host == phiphi.t2.* ]]; then
  dir=/data/nfsmic/scratch/toymc
  n_sim_thr=12
else
  dir=/tmp/${USER}/toymc
  n_sim_thr=8
fi
