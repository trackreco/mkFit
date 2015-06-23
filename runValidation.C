void runValidation(){
  gROOT->LoadMacro("PlotValidation.cpp++g");

  // First is input name of root file
  // Second is output name of directory/rootfile/file plots
  // Third is output type of plots
  
  PlotValidation Validation("valtree.root","seedval_conf_inward_ENDTOEND_adjustedIII","pdf");
  Validation.Validation(Bool_t(true)); // bool to move input root file to output directory
}
