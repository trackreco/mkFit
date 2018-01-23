#! /bin/bash

# in case this is run alone
source xeon_scripts/common_variables.sh

dir=${1:-benchmarks}

# Move benchmark plots: subroutine build benchmarks
builddir="Benchmarks"
mkdir -p ${dir}/${builddir}
mkdir -p ${dir}/${builddir}/logx

for arch in SNB KNL
do
    for benchmark in TH VU
    do
	oBase=${arch}_${physics_sample}_${benchmark}

	mv ${oBase}_"time".png ${dir}/${builddir}
	mv ${oBase}_"speedup".png ${dir}/${builddir}

	mv ${oBase}_"time_logx".png ${dir}/${builddir}/logx
	mv ${oBase}_"speedup_logx".png ${dir}/${builddir}/logx
    done
done

# Move multiple events in flight plots
meifdir="MultEvInFlight"
mkdir -p ${dir}/${meifdir}
mkdir -p ${dir}/${meifdir}/logx

for arch in SNB KNL
do
    for build in CE FV ; do
        oBase=${arch}_${physics_sample}_${build}_"MEIF"

        mv ${oBase}_"time".png ${dir}/${meifdir}
        mv ${oBase}_"speedup".png ${dir}/${meifdir}

        mv ${oBase}_"time_logx".png ${dir}/${meifdir}/logx
        mv ${oBase}_"speedup_logx".png ${dir}/${meifdir}/logx
    done
done

# Move plots from text dump
dumpdir="PlotsFromDump"
mkdir -p ${dir}/${dumpdir}
mkdir -p ${dir}/${dumpdir}/diffs

for build in BH STD CE FV
do
    for var in nHits pt eta phi
    do
	mv ${physics_sample}_${build}_${var}.png ${dir}/${dumpdir}
	mv ${physics_sample}_${build}_"d"${var}.png ${dir}/${dumpdir}/diffs
    done
done

# Move ROOT validation
rootdir="ROOTVAL"
mkdir -p ${dir}/${rootdir}

for build in BH STD CE FV
do
    vBase=${val_arch}_${rootval_sample}_${build}
    mv validation_${vBase}_"ROOTVAL"/totals_validation_${vBase}_"ROOTVAL".txt ${dir}/${rootdir}
done

for rate in eff ineff_barrel ineff_endcap dr fr 
do
    for var in pt phi eta
    do 
	for pt in 0.0 0.9 2.0
	do
	    mv ${val_arch}_${rootval_sample}_${rate}_${var}_"build"_pt${pt}_"ROOTVAL".png ${dir}/${rootdir}
	done
    done
done

# Move CMSSW validation
cmsswdir="CMSSWVAL"
mkdir -p ${dir}/${cmsswdir}

for build in BH STD CE FV
do
    vBase=${val_arch}_${physics_sample}_${build}
    mv validation_${vBase}_"CMSSWVAL"/totals_validation_${vBase}_"CMSSWVAL"_cmssw.txt ${dir}/${cmsswdir}
done

for trk in build fit
do
    mkdir -p ${dir}/${cmsswdir}/${trk}
    mkdir -p ${dir}/${cmsswdir}/${trk}/diffs
done

for rate in eff ineff_barrel ineff_endcap dr fr 
do
    for var in pt phi eta
    do
	for trk in build fit
	do
	    for pt in 0.0 0.9 2.0
	    do
		mv ${val_arch}_${physics_sample}_${rate}_${var}_${trk}_pt${pt}_"CMSSWVAL".png ${dir}/${cmsswdir}/${trk}
	    done
	done
    done
done    

for coll in bestmatch allmatch
do 
    for var in nHits invpt phi eta
    do
	for trk in build fit
	do
	    for pt in 0.0 0.9 2.0
	    do
		mv ${val_arch}_${physics_sample}_${coll}_d${var}_${trk}_pt${pt}_"CMSSWVAL".png ${dir}/${cmsswdir}/${trk}/diffs
	    done
	done
    done
done
