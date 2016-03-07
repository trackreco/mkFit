#! /bin/bash

[ -z "$ROOTSYS" ] && source ~matevz/root/bin/thisroot.sh 

source benchmark.sh
source benchmark-fit.sh  

python makeBenchmarkPlots.py host

root -b -q -l makeBenchmarkPlots.C
root -b -q -l makeBenchmarkPlotsFit.C

source benchmark-mic.sh  
source benchmark-mic-fit.sh  

python makeBenchmarkPlots.py mic

sed -i 's/bool isMic = false/bool isMic = true/g' makeBenchmarkPlots.C
sed -i 's/bool isMic = false/bool isMic = true/g' makeBenchmarkPlotsFit.C
root -b -q -l makeBenchmarkPlots.C
root -b -q -l makeBenchmarkPlotsFit.C
sed -i 's/bool isMic = true/bool isMic = false/g' makeBenchmarkPlots.C
sed -i 's/bool isMic = true/bool isMic = false/g' makeBenchmarkPlotsFit.C

python makePlotsFromDump.py _host_10x20k_BH_NVU1_NTH1
python makePlotsFromDump.py _host_10x20k_BH_NVU8int_NTH21
python makePlotsFromDump.py _mic_10x20k_BH_NVU1_NTH1
python makePlotsFromDump.py _mic_10x20k_BH_NVU16int_NTH210
sed -i 's/\/\/TString test = "BH"/TString test = \"BH\"/g' makePlotsFromDump.C
root -b -q -l makePlotsFromDump.C
sed -i 's/TString test = "BH"/\/\/TString test = \"BH\"/g' makePlotsFromDump.C

python makePlotsFromDump.py _host_10x20k_CE_NVU1_NTH1
python makePlotsFromDump.py _host_10x20k_CE_NVU8int_NTH21
python makePlotsFromDump.py _mic_10x20k_CE_NVU1_NTH1
python makePlotsFromDump.py _mic_10x20k_CE_NVU16int_NTH210
sed -i 's/\/\/TString test = "CE"/TString test = \"CE\"/g' makePlotsFromDump.C
root -b -q -l makePlotsFromDump.C
sed -i 's/TString test = "CE"/\/\/TString test = \"CE\"/g' makePlotsFromDump.C

python makePlotsFromDump.py _host_10x20k_CEST_NVU1_NTH1
python makePlotsFromDump.py _host_10x20k_CEST_NVU8int_NTH21
python makePlotsFromDump.py _mic_10x20k_CEST_NVU1_NTH1
python makePlotsFromDump.py _mic_10x20k_CEST_NVU16int_NTH210
sed -i 's/\/\/TString test = "CEST"/TString test = \"CEST\"/g' makePlotsFromDump.C
root -b -q -l makePlotsFromDump.C
sed -i 's/TString test = "CEST"/\/\/TString test = \"CEST\"/g' makePlotsFromDump.C

python makePlotsFromDump.py _host_10x20k_ST_NVU1_NTH1
python makePlotsFromDump.py _host_10x20k_ST_NVU8int_NTH21
python makePlotsFromDump.py _mic_10x20k_ST_NVU1_NTH1
python makePlotsFromDump.py _mic_10x20k_ST_NVU16int_NTH210
sed -i 's/\/\/TString test = "ST"/TString test = \"ST\"/g' makePlotsFromDump.C
root -b -q -l makePlotsFromDump.C
sed -i 's/TString test = "ST"/\/\/TString test = \"ST\"/g' makePlotsFromDump.C
