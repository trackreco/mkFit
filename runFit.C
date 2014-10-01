void runFit(){
  gROOT->LoadMacro("PlotFit.cpp++g");

  // First is input name of root file
  // Second is output name of directory/rootfile/file plots
  // Third is output type of plots
  
  // Can make two objects to do fitting plots at same time for two different releases
  // Can run comparator immediately following running both objects

  PlotFit fitPlotOld("validationtree_cylinder.root","cylinderFit","pdf");
  fitPlotOld.PlotAllHistos();

  PlotFit fitPlotNew("validationtree_poly.root","polyFit","pdf");
  fitPlotNew.PlotAllHistos();
}
