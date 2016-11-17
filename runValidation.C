void setupcpp11()
{ // customize ACLiC's behavior ...
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

void runValidation(TString test = "", Bool_t mvInput = true, Bool_t fullval = false,
		   Bool_t saveAs = false, TString image = "pdf")
{
  setupcpp11(); //  use this to get PlotValidation to compile ... phiphi ROOT build has ACLiC with C++98!

  gROOT->LoadMacro("PlotValidation.cpp++g");

  // PlotValidation arguments
  // First is additional input name of root file
  // Second is output name of directory/rootfile/file plots
  // The first boolean argument == true to move input root file to output directory, false to keep input file where it is.
  // Sedcond bool == true for full validation, false for only "main" validation.
  // Main validation includes efficiencies, fake rates, duplicate rates. Also momentum pulls, nHits plots, timing plots.
  // Full validation includes main validation plus geometry plots, positions pulls, and CF pulls/residuals.  
  // Full validation also includes branching plots and segmenting plots
  // Third Bool is saving the image files
  // Last argument is output type of plots


  PlotValidation Val(Form("valtree%s.root",test.Data()),Form("validation%s",test.Data()),
		     mvInput,fullval,saveAs,image);
  Val.Validation(); 
}
