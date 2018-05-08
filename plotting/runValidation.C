#include "PlotValidation.hh"
#include "PlotValidation.cpp"

void runValidation(TString test = "", Bool_t computePulls = false, Bool_t cmsswComp = false, 
		   Bool_t mvInput = true, Bool_t saveAs = false, TString image = "pdf")
{
  // PlotValidation arguments
  // First is additional input name of root file
  // Second is name of output directory
  // First boolean argument is to compute momentum pulls: currently implemented only when sim track states are available!
  // Second boolean argument is to do special CMSSW validation
  // The third boolean argument == true to move input root file to output directory, false to keep input file where it is.
  // Fourth Bool is saving the image files
  // Last argument is output type of plots

  PlotValidation Val(Form("valtree%s.root",test.Data()),Form("validation%s",test.Data()),computePulls,cmsswComp,mvInput,saveAs,image);
  Val.Validation(); 
}
