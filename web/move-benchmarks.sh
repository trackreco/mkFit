#!/bin/bash

dir=${1:-benchmarks}
base=CMSSW_TTbar_PU70

echo "Moving plots and text files locally to ${outdir}"

# Move benchmark plots: subroutine build benchmarks
builddir="Benchmarks"
mkdir -p ${dir}/${builddir}
mkdir -p ${dir}/${builddir}/logx

for arch in SNB KNC KNL
do
    for benchmark in TH VU
    do
	oBase=${arch}_${base}_${benchmark}

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

for arch in SNB KNC KNL
do
    oBase=${arch}_${base}_"MEIF"

    mv ${oBase}_"time".png ${dir}/${meifdir}
    mv ${oBase}_"speedup".png ${dir}/${meifdir}

    mv ${oBase}_"time_logx".png ${dir}/${meifdir}/logx
    mv ${oBase}_"speedup_logx".png ${dir}/${meifdir}/logx
done

# Move plots from text dump
dumpdir="PlotsFromDump"
mkdir -p ${dir}/${dumpdir}

for build in BH STD CE
do
    mv ${base}_${build}_*"nHits".png ${dir}/${dumpdir}
    mv ${base}_${build}_*"pt".png ${dir}/${dumpdir}
    mv ${base}_${build}_*"eta".png ${dir}/${dumpdir}
    mv ${base}_${build}_*"phi".png ${dir}/${dumpdir}
done

# Move ROOT validation
rootdir="ROOTVAL"
mkdir -p ${dir}/${rootdir}

for build in BH STD CE
do
    vBase=SNB_${base}_${build}
    mv validation_${vBase}_"ROOTVAL"/totals_validation_${vBase}_"ROOTVAL".txt ${dir}/${rootdir}
done
mv SNB_${base}_*_"ROOTVAL".png ${dir}/${rootdir}

# Move CMSSW validation
cmsswdir="CMSSWVAL"
mkdir -p ${dir}/${cmsswdir}

for build in BH STD CE
do
    vBase=SNB_${base}_${build}
    mv validation_${vBase}_"CMSSWVAL"/totals_validation_${vBase}_"CMSSWVAL"_cmssw.txt ${dir}/${cmsswdir}
done
mv SNB_${base}_*_"CMSSWVAL".png ${dir}/${cmsswdir}

host=kmcdermo@lxplus.cern.ch
whost=${host}":~/www"
echo "Moving plots and text files remotely to ${whost}"
scp -r ${dir} ${whost}

echo "Executing remotely ./makereadable.sh ${outdir}"
ssh ${host} bash -c "'
cd www
./makereadable.sh ${outdir}
exit
'"

echo "Removing local files"
rm -rf log_*.txt
rm -rf *.root
rm -rf *.png
rm -rf validation_*

rm -rf ${dir}
