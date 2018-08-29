#include "plotting/PlotsFromDump.cpp+"

void makePlotsFromDump(const TString & sample, const TString & build, const TString & suite)
{
  PlotsFromDump Plots(sample,build,suite);
  Plots.RunPlotsFromDump();
}
