#include "PlotMEIFBenchmarks.hh"
#include "PlotMEIFBenchmarks.cpp"

void makeMEIFBenchmarkPlots(const TString & arch, const TString & sample, const TString & build)
{
  PlotMEIFBenchmarks MEIFBenchmarks(arch,sample,build);
  MEIFBenchmarks.RunMEIFBenchmarkPlots();
}
