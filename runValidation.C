void runValidation(){
  gROOT->LoadMacro("PlotValidation.cpp++g");

  // First is input name of root file
  // Second is output name of directory/rootfile/file plots
  // Third is output type of plots
  
  PlotValidation Validation("valtree.root","mcENDTOEND","pdf");
  Validation.Validation(Bool_t(false));
}
