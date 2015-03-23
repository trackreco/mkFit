void runSeedEff(){
  gROOT->LoadMacro("PlotSeedEff.cpp++g");

  // First is input name of root file
  // Second is output name of directory/rootfile/file plots
  // Third is output type of plots

  PlotSeedEff seed500("valtree.root","seed500const","pdf");
  seed500.PlotAllHistos();

  //  PlotSeedEff seed2000("valtree_seed2000.root","seed2000","pdf");
  // seed2000.PlotAllHistos();

  // PlotSeedEff seed5000("valtree_seed5000.root","seed5000","pdf");
  // seed5000.PlotAllHistos();
}
