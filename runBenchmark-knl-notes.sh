./benchmark-knl.sh
./benchmark-knl-fit.sh

# Run where root exists

python makeBenchmarkPlots.py knl

root -b -q -l makeBenchmarkPlots.C\(0,0,0,1\)
root -b -q -l makeBenchmarkPlotsFit.C\(0,0,1\)

# Run parts of this where combined with mic/host.
# Need to decide what we compare against ... or it will be a mess.
for test in BH CE CEST ST TBBST; do

    # Orig
    python makePlotsFromDump.py _host_10x20k_${test}_NVU1_NTH1
    python makePlotsFromDump.py _host_10x20k_${test}_NVU8int_NTH21
    python makePlotsFromDump.py _mic_10x20k_${test}_NVU1_NTH1
    python makePlotsFromDump.py _mic_10x20k_${test}_NVU16int_NTH210

    # KNL
    python makePlotsFromDump.py _knl_10x20k_${test}_NVU1_NTH1
    python makePlotsFromDump.py _knl_10x20k_${test}_NVU16int_NTH210

    root -b -q -l makePlotsFromDump.C\(\"${test}\"\)
done
