#include "PlotBenchmarks.hh"
#include "PlotBenchmarks.cpp"

void makeBenchmarkPlots(const TString & arch, const TString & sample)
{
  PlotBenchmarks Benchmarks(arch,sample);
  Benchmarks.RunBenchmarkPlots();
}
