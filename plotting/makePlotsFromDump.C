#if  !defined(__CINT__)
#include "plotting/PlotsFromDump.hh"
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

void makePlotsFromDump(const TString & sample, const TString & build)
{
  setupcpp11(); //  use this to get PlotsFromDump to compile ... phiphi ROOT build has ACLiC with C++98!

  gROOT->LoadMacro("plotting/PlotsFromDump.cpp+g");

  PlotsFromDump Plots(sample,build);
  Plots.RunPlotsFromDump();
}
