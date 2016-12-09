#! /bin/bash

[ -z "$ROOTSYS" ] && source ~matevz/root/bin/thisroot.sh

##### Generate toyMC samples first #####
./generateToyMCsamples.sh

##### SNB Tests #####
## ToyMC ##
./benchmark-snb-toymc-barrel-build.sh
./benchmark-snb-toymc-barrel-fit.sh
python makeBenchmarkPlots.py SNB ToyMC Barrel
root -b -q -l makeBenchmarkPlots.C
root -b -q -l makeBenchmarkPlotsFit.C

./benchmark-snb-toymc-endcap-fit.sh
python makeBenchmarkPlots.py SNB ToyMC Endcap
root -b -q -l makeBenchmarkPlotsFit.C\(0,1\)

## CMSSW ##
./benchmark-snb-cmssw-barrel-build.sh
python makeBenchmarkPlots.py SNB CMSSW Barrel
root -b -q -l makeBenchmarkPlots.C\(0,1\)

./benchmark-snb-cmssw-endcap-build.sh
python makeBenchmarkPlots.py SNB CMSSW Endcap
root -b -q -l makeBenchmarkPlots.C\(0,1,1\)

##### KNC Tests #####
## ToyMC ##
./benchmark-knc-toymc-barrel-build.sh
./benchmark-knc-toymc-barrel-fit.sh
python makeBenchmarkPlots.py KNC ToyMC Barrel
root -b -q -l makeBenchmarkPlots.C\(1\)
root -b -q -l makeBenchmarkPlotsFit.C\(1\)

./benchmark-knc-toymc-endcap-fit.sh
python makeBenchmarkPlots.py KNC ToyMC Endcap
root -b -q -l makeBenchmarkPlotsFit.C\(1,1\)

## CMSSW ##
./benchmark-knc-cmssw-barrel-build.sh
python makeBenchmarkPlots.py KNC CMSSW Barrel
root -b -q -l makeBenchmarkPlots.C\(1,1\)

##### nHits plots #####
for test in BH STD CE; do 
    echo "Making nHits plots for ToyMC Barrel:" ${test}
    python makePlotsFromDump.py SNB ToyMC Barrel ${test} NVU1_NTH1
    python makePlotsFromDump.py SNB ToyMC Barrel ${test} NVU8int_NTH24
    python makePlotsFromDump.py KNC ToyMC Barrel ${test} NVU1_NTH1
    python makePlotsFromDump.py KNC ToyMC Barrel ${test} NVU16int_NTH240
    root -b -q -l makePlotsFromDump.C\(\"${test}\"\)
    
    echo "Making nHits plots for CMSSW Barrel:" ${test}
    python makePlotsFromDump.py SNB CMSSW Barrel ${test} NVU1_NTH1
    python makePlotsFromDump.py SNB CMSSW Barrel ${test} NVU8int_NTH24
    python makePlotsFromDump.py KNC CMSSW Barrel ${test} NVU1_NTH1
    python makePlotsFromDump.py KNC CMSSW Barrel ${test} NVU16int_NTH240
    root -b -q -l makePlotsFromDump.C\(\"${test}\",1\)
     
    echo "Making nHits plots for CMSSW Endcap:" ${test}
    python makePlotsFromDump.py SNB CMSSW Endcap ${test} NVU1_NTH1
    python makePlotsFromDump.py SNB CMSSW Endcap ${test} NVU8int_NTH24
    root -b -q -l makePlotsFromDump.C\(\"${test}\",1,1\)
done

##### Validation tests #####
./validation-snb-toymc-barrel-build.sh
for test in BH STD CE; do
    root -b -q -l runValidation.C\(\"_SNB_ToyMC_Barrel_${test}\"\)
done
root -b -q -l makeValidation.C\(\"SNB_ToyMC_Barrel\"\)

make distclean
