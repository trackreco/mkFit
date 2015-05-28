void runValidation(){
  gROOT->LoadMacro("PlotValidation.cpp++g");

  // First is input name of root file
  // Second is output name of directory/rootfile/file plots
  // Third is output type of plots
  
  PlotValidation Validation("valtree_realseeds.root","seedval","pdf");
  Validation.Validation(Bool_t(true));
  //Validation.PlotNumerNHits();
}
