#if  !defined(__CINT__)
#include "plotting/StackValidation.hh"
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

void makeValidation(const TString & label = "", const TString & extra = "", const Bool_t cmsswComp = false)
{
  setupcpp11(); //  use this to get StackValidation to compile ... phiphi ROOT build has ACLiC with C++98!

  gROOT->LoadMacro("plotting/StackValidation.cpp+g");

  StackValidation Stacks(label,extra,cmsswComp);
  Stacks.MakeValidationStacks();
}
