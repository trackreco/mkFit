#include "plotting/PlotsFromDump.cpp+"

void makePlotsFromDump(const TString & sample, const TString & build)
{
  PlotsFromDump Plots(sample,build);
  Plots.RunPlotsFromDump();
}
