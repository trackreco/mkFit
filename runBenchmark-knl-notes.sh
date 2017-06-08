./benchmark-knl-toymc-barrel-build.sh
./benchmark-knl-toymc-barrel-fit.sh

# Run where root exists

python makeBenchmarkPlots.py KNL

root -b -q -l makeBenchmarkPlots.C\(0,0,0,1\)
root -b -q -l makeBenchmarkPlotsFit.C\(0,0,1\)

# Run parts of this where combined with mic/host.
# Need to decide what we compare against ... or it will be a mess.
for test in BH STD CE; do

    python makePlotsFromDump.py KNL ToyMC Barrel ${test} NVU1_NTH1
    python makePlotsFromDump.py KNL ToyMC Barrel ${test} NVU16int_NTH210

    root -b -q -l makePlotsFromDump.C\(\"${test}\"\)
done
