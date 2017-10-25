#! /bin/bash

# To be sourced where needed

host=`hostname`
user=${TOYMCUSER:-$USER}

if [[ $host == phi2.t2.* ]]; then
  dir=/data/scratch/$user/toymc
  n_sim_thr=128
elif [[ $host == phiphi.t2.* ]]; then
  dir=/data/nfsmic/$user/toymc
  n_sim_thr=12
else
  dir=/tmp/$user/toymc
  n_sim_thr=8
fi
