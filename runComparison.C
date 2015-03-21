void runComparison(){
  gROOT->LoadMacro("PlotCompare.cpp++g");

  // Arguments of PlotComparator constructor 
  // 1+2 are names of strings to be used to label old/new releases for comparison
  // Third argument output string for directory name/file names 
  // Fourth argument for output type for plots

  // First argument of PlotCompare are list of plots in each root file to be compared
  // Last two arguments are those root files, first is old, second is new release 
  //  plotter.PlotCompare("fit.list","XYZ_SCATTER_fit/XYZ_SCATTER_fit.root","SOLID_SMEAR_fit/SOLID_SMEAR_fit.root");
  //  plotter.PlotCompare("build.list","build_validationtree_XYZ_SCATTER.root","build_validationtree_5000.root");
    
  PlotComparator plotter5000("phiSegOnly5000","phiEtaSeg5000","phiSegOnly_vs_phiEtaSeg_5000","pdf");
  plotter5000.PlotCompare("buildEff.list","phionly/5000phionly/5000phionly.root","phieta/5000phieta/5000phieta.root");
 
  PlotComparator plotter500("phiSegOnly500","phiEtaSeg500","phiSegOnly_vs_phiEtaSeg_500","pdf");
  plotter500.PlotCompare("buildEff.list","phionly/500phionly/500phionly.root","phieta/500phieta/500phieta.root");
 
  PlotComparator plotter2000("phiSegOnly2000","phiEtaSeg2000","phiSegOnly_vs_phiEtaSeg_2000","pdf");
  plotter2000.PlotCompare("buildEff.list","phionly/2000phionly/2000phionly.root","phieta/2000phieta/2000phieta.root");
  /*  
  PlotComparator plotter5001("500","1","500_vs_1_trksPerEvt","pdf");
  plotter5001.PlotCompare("buildEff.list","500_buildEff/500_buildEff.root","1_buildEff/1_buildEff.root");
  
  PlotComparator plotter50001("5000","1","5000_vs_1_trksPerEvt","pdf");
  plotter50001.PlotCompare("buildEff.list","5000_buildEff/5000_buildEff.root","1_buildEff/1_buildEff.root");*/
}
