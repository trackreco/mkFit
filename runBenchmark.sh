#! /bin/bash

[ -z "$ROOTSYS" ] && source ~matevz/root/bin/thisroot.sh

##### Generate toyMC samples first #####
./generateToyMCsamples.sh

##### SNB Tests #####
./benchmark-snb-toymc-barrel-build.sh
./benchmark-snb-toymc-barrel-fit.sh
python makeBenchmarkPlots.py snb
root -b -q -l makeBenchmarkPlots.C
root -b -q -l makeBenchmarkPlotsFit.C

./benchmark-snb-toymc-endcap-fit.sh
python makeBenchmarkPlots.py snb_endcap
root -b -q -l makeBenchmarkPlotsFit.C\(0,1\)

./benchmark-snb-cmssw-barrel-build.sh
python makeBenchmarkPlots.py snb cmssw
root -b -q -l makeBenchmarkPlots.C\(0,1\)

##### KNC Tests #####
./benchmark-knc-toymc-barrel-build.sh
./benchmark-knc-toymc-barrel-fit.sh
python makeBenchmarkPlots.py knc
root -b -q -l makeBenchmarkPlots.C\(1\)
root -b -q -l makeBenchmarkPlotsFit.C\(1\)

./benchmark-knc-toymc-endcap-fit.sh
python makeBenchmarkPlots.py knc_endcap
root -b -q -l makeBenchmarkPlotsFit.C\(1,1\)

./benchmark-knc-cmssw-barrel-build.sh
python makeBenchmarkPlots.py knc cmssw
root -b -q -l makeBenchmarkPlots.C\(1,1\)

##### nHits plots #####
for test in BH CE; do # add in STD once standard building converted to TBB
    echo "Making nHits plots for ToyMC barrel:" ${test}
    python makePlotsFromDump.py _snb_20x10k_${test}_NVU1_NTH1
    python makePlotsFromDump.py _snb_20x10k_${test}_NVU8int_NTH24
    python makePlotsFromDump.py _knc_20x10k_${test}_NVU1_NTH1
    python makePlotsFromDump.py _knc_20x10k_${test}_NVU16int_NTH240
    root -b -q -l makePlotsFromDump.C\(\"${test}\"\)

    echo "Making nHits plots for CMSSW barrel:" ${test}
    python makePlotsFromDump.py _snb_100xTTbarPU35_${test}_NVU1_NTH1
    python makePlotsFromDump.py _snb_100xTTbarPU35_${test}_NVU8int_NTH24
    python makePlotsFromDump.py _knc_100xTTbarPU35_${test}_NVU1_NTH1
    python makePlotsFromDump.py _knc_100xTTbarPU35_${test}_NVU16int_NTH240
    root -b -q -l makePlotsFromDump.C\(\"${test}\",1\)
done

##### Validation tests #####
./validation-snb-toymc-barrel-build.sh
root -b -q -l runValidation.C\(\"_BH\"\)
#root -b -q -l runValidation.C\(\"_STD\"\)
root -b -q -l runValidation.C\(\"_CE\"\)
root -b -q -l makeValidation.C

make distclean
