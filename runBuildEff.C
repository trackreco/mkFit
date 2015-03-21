void runBuildEff(){
  gROOT->LoadMacro("PlotBuildEff.cpp++g");

  // First is input name of root file
  // Second is output name of directory/rootfile/file plots
  // Third is output type of plots
  // Fourth is boolean to move input to same output directory

  PlotBuildEff buildEffPlots("valtree.root","ss_PE_ED2.0_BS_r0.0025_z1.0_1tk_scatter0.0025","pdf",true);
  buildEffPlots.Validation();

}
