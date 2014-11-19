void runBuildEff(){
  gROOT->LoadMacro("PlotBuildEff.cpp++g");

  // First is input name of root file
  // Second is output name of directory/rootfile/file plots
  // Third is output type of plots
  
  PlotBuildEff buildPlots_Old("build_validationtree.root","function","pdf");
  buildPlots_Old.PlotAllHistos();
  /*
  PlotBuildEff buildPlots_New("build_validationtree5000.root","5000_buildEff","pdf");
  buildPlots_New.PlotAllHistos();

  PlotBuildEff buildPlots_New1("build_validationtree1.root","1_buildEff","pdf");
  buildPlots_New1.PlotAllHistos();
  */
}
