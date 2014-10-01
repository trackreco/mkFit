#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

#include <iostream>
#include <fstream>
#include <cmath>

class PlotComparator
{
 public:
  
  PlotComparator(TString oldName, TString newName, TString outName, TString outType);
  ~PlotComparator();
  void SetUpPlotter();
  void CreatePlot(TH1F * hist_old, TH1F * hist_new);
  void PlotCompareFit(TString listOfPlots, TString inRootOldName, TString inRootNewName);
  void PlotCompareBuild(TString listOfPlots, TString inRootOldName, TString inRootNewName); 
  void PlotLister(TString listOfPlots);

private:
  TCanvas * fCompareCanvas;
  TFile * fRootOut;

  UInt_t fNPlots;
  TString * fPlotsList;

  TString fOldName;
  TString fNewName;
  TString fOutName;
  TString fOutType;
};
