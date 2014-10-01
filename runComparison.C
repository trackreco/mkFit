void runComparison(){
  gROOT->LoadMacro("PlotCompare.cpp++g");

  // Arguments 1+2 are names of strings to be used to label old/new releases for comparison
  // Third argument output string for directory name/file names 
  // Fourth argument for output type for plots
  
  PlotComparator plotter("cyl","poly","CylVsPoly","pdf");

  // Two functions to compare fitter and builder
  // First argument are list of plots in each root file to be compared
  // Last two arguments are those root files, first is old, second is new release 
  plotter.PlotCompareFit("fit.list","cylinderFit/cylinderFit.root","polyFit/polyFit.root");
  plotter.PlotCompareBuild("build.list","build_validationtree_cylinder.root","build_validationtree_poly.root");
}
