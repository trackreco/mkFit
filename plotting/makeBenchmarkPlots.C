#if  !defined(__CINT__)
#include "plotting/PlotBenchmarks.hh"
#endif

void setupcpp11()
{ 
  // customize ACLiC's behavior ...
  TString o;
  // customize MakeSharedLib
  o = TString(gSystem->GetMakeSharedLib());
  o = o.ReplaceAll(" -c ", " -std=c++0x -c ");
  gSystem->SetMakeSharedLib(o.Data());
  // customize MakeExe
  o = TString(gSystem->GetMakeExe());
  o = o.ReplaceAll(" -c ", " -std=c++0x -c ");
  gSystem->SetMakeExe(o.Data());
} 

void makeBenchmarkPlots(const TString & arch, const TString & sample)
{
  setupcpp11(); //  use this to get PlotBenchmarks to compile ... phiphi ROOT build has ACLiC with C++98!

  gROOT->LoadMacro("plotting/PlotBenchmarks.cpp+g");

  PlotBenchmarks Benchmarks(arch,sample);
  Benchmarks.RunBenchmarkPlots();
}
