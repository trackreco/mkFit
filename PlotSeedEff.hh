#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TPaveText.h"

#include <iostream>
#include <cmath>

class PlotSeedEff
{
public:
  PlotSeedEff(TString inName, TString outName, TString outType);
  ~PlotSeedEff();
  void PlotAllHistos();
  void PlotEffandExtra();
  void SetUpRatioPlots(TH1F *& plot, TString name, TString title, Int_t nBins, Double_t lowEdge, Double_t upEdge, TString xtitle);
  void CalculateRatioPlot(TH1F * numer, TH1F * denom, TH1F *& ratioPlot, UInt_t type);
  //  void PlaceNEvents(TH1F * numer, TH1F * denom, TH1F *& ratioPlot);
  void DrawSaveTH1Plot(TH1F * histo, TString outputName);

private:
  TString fOutType;
  TString fOutName;
  TFile * fOutRoot;
  TFile * fInRoot;
  TCanvas * fTH1Canv;
  //  TPaveText * fTPave; 
};
