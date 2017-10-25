#! /bin/bash

make -j 12

. data-dir-location.sh

mkdir -p ${dir}

# Building test: Barrel
if [ ! -f ${dir}/simtracks_barrel_20x10k.bin ]; then
    echo "++++Generating 10k tracks/event * 20 events for ToyMC barrel-only building tests++++"
    ./mkFit/mkFit --num-thr-sim ${n_sim_thr} --num-events 20 --num-tracks 10000 --write --file-name simtracks_barrel_20x10k.bin
    mv simtracks_barrel_20x10k.bin ${dir}/
fi

# Fitting test: Barrel
if [ ! -f ${dir}/simtracks_barrel_1kx10k.bin ]; then
    echo "++++Generating 10k tracks/event * 1k events for ToyMC barrel-only fiting tests++++"
    ./mkFit/mkFit --num-thr-sim ${n_sim_thr} --num-events 1000 --num-tracks 10000 --write --file-name simtracks_barrel_1kx10k.bin
    mv simtracks_barrel_1kx10k.bin ${dir}/
fi

# Fitting test: Endcap
if [ ! -f ${dir}/simtracks_endcap_1kx10k.bin ]; then
    echo "++++Generating 10k tracks/event * 1k events for ToyMC endcap-only fiting tests++++"
    ./mkFit/mkFit --num-thr-sim ${n_sim_thr} --endcap-test --num-events 1000 --num-tracks 10000 --write --file-name simtracks_endcap_1kx10k.bin
    mv simtracks_endcap_1kx10k.bin ${dir}/
fi

make clean

# Validation tests: Barrel
if [ ! -f ${dir}/simtracks_fulldet_400x2p5k_val.bin ]; then
    echo "++++Generating 2.5k tracks/event * 400 events for ToyMC barrel-only validation tests++++"
    sed -i 's/#WITH_ROOT := yes/WITH_ROOT := yes/g' Makefile.config
    make -j 12
    ./mkFit/mkFit --num-thr-sim 12 --root-val --num-events 400 --num-tracks 2500 --write --file-name simtracks_fulldet_400x2p5k_val.bin
    mv simtracks_fulldet_400x2p5k_val.bin ${dir}/
    sed -i 's/WITH_ROOT := yes/#WITH_ROOT := yes/g' Makefile.config
fi

make clean
