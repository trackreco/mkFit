#!/bin/bash

vhb_v=(10)
mhp_v=(10)
mhc_v=(3 4)
mch_v=(1 2 99)
c2c_v=(15 30 45 60)

for vhb in "${vhb_v[@]}"
do
    echo ${vhb};
    sed -i "s/validHitBonus_ = XX/validHitBonus_ = ${vhb}/g" Config.h
    for mhp in "${mhp_v[@]}"
    do
	echo ${mhp};
	sed -i "s/missingHitPenalty_ = XX/missingHitPenalty_ = ${mhp}/g" Config.h
	for mhc in "${mhc_v[@]}"
	do
	    echo ${mhc};
	    sed -i "s/maxHolesPerCand  = XX/maxHolesPerCand  = ${mhc}/g" Geoms/CMS-2017.cc

	    for mch in "${mch_v[@]}"
	    do
		echo ${mch};
		sed -i "s/maxConsecHoles   = XX/maxConsecHoles   = ${mch}/g" Geoms/CMS-2017.cc

		for c2c in "${c2c_v[@]}"
		do
		    echo ${c2c};
		    sed -i "s/chi2Cut          = XX/chi2Cut          = ${c2c}/g" Geoms/CMS-2017.cc
		    grep validHitBonus Config.h
		    grep missingHitPenalty Config.h
		    grep maxHolesPerCand Geoms/CMS-2017.cc
		    grep maxConsecHoles  Geoms/CMS-2017.cc
		    grep chi2Cut Geoms/CMS-2017.cc
		    make -j 32 WITH_ROOT:=1 AVX_512:=1
		    ./mkFit/mkFit --silent --cmssw-n2seeds --num-thr 64 --num-thr-ev 32 --input-file /data2//pu50-ccc-hs.bin --num-events 500 --remove-dup --sim-val --try-to-save-sim-info --backward-fit-pca --mtv-require-seeds --build-ce
		    hadd -O valtree-SKL-SP_CMSSW_TTbar_PU50_CE_SIMVALSEED-test.root valtree_*.root
		    rm valtree_*.root
		    mv valtree-SKL-SP_CMSSW_TTbar_PU50_CE_SIMVALSEED-test.root valtree-SKL-SP_CMSSW_TTbar_PU50_CE_SIMVALSEED-test_${vhb}_${mhp}_${mhc}_${mch}_${c2c}.root
		    sed -i "s/chi2Cut          = ${c2c}/chi2Cut          = XX/g" Geoms/CMS-2017.cc
		done
		sed -i "s/maxConsecHoles   = ${mch}/maxConsecHoles   = XX/g" Geoms/CMS-2017.cc
	    done
	    sed -i "s/maxHolesPerCand  = ${mhc}/maxHolesPerCand  = XX/g" Geoms/CMS-2017.cc
	done
	sed -i "s/missingHitPenalty_ = ${mhp}/missingHitPenalty_ = XX/g" Config.h
    done
    sed -i "s/validHitBonus_ = ${vhb}/validHitBonus_ = XX/g" Config.h
done
