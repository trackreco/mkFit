void setupcpp11(){ // customize ACLiC's behavior ...
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

void runValidation() {
  setupcpp11(); //  use this to get PlotValidation to compile ... phiphi ROOT build has ACLiC with C++98!

  gROOT->LoadMacro("PlotValidation.cpp++g");

  // First is input name of root file
  // Second is output name of directory/rootfile/file plots
  // Third is output type of plots

  // Validation() has two boolean arguments
  // First bool == true for full validation, false for only "main" validation.
  // Main validation includes efficiencies, fake rates, duplicate rates. Also momentum pulls, nHits plots, timing plots.
  // Full validation includes main validation plus geometry plots, positions pulls, and CF pulls/residuals.  
  // Full validation also includes branching plots and segmenting plots
  // The second boolean argument == true to move input root file to output directory, false to keep input file where it is.

  PlotValidation Val("valtree.root","full_validation","pdf");
  Val.Validation( Bool_t(false), Bool_t(true) ); 
}
