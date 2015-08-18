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
  
  PlotValidation Validation("valtree.root","5000tk_100evts_630phi_10eta","pdf");
  Validation.Validation(Bool_t(true)); // bool to move input root file to output directory
}
