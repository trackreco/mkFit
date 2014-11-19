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
    
  PlotComparator plotter500f("500","f","500_vs_f","pdf");
  plotter500f.PlotCompare("buildEff.list","500_buildEff/500_buildEff.root","function/function.root");
  /*  
  PlotComparator plotter5001("500","1","500_vs_1_trksPerEvt","pdf");
  plotter5001.PlotCompare("buildEff.list","500_buildEff/500_buildEff.root","1_buildEff/1_buildEff.root");
  
  PlotComparator plotter50001("5000","1","5000_vs_1_trksPerEvt","pdf");
  plotter50001.PlotCompare("buildEff.list","5000_buildEff/5000_buildEff.root","1_buildEff/1_buildEff.root");*/
}
