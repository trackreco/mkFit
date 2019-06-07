#include "plotting/PlotsFromDump.cpp+"

void makePlotsFromDump(const TString & sample, const TString & build, const TString & suite, const int useLNX)
{
  PlotsFromDump Plots(sample,build,suite,useLNX);
  Plots.RunPlotsFromDump();
}
