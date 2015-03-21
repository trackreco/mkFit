#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

#include <iostream>
#include <fstream>
#include <cmath>

class PlotBuildEff
{
public:
  PlotBuildEff(TString inName, TString outName, TString outType, Bool_t mvIn);
  ~PlotBuildEff();
  void Validation();
  void PlotEffandFakeRate();
  void SetUpRatioPlots(TH1F *& plot, TString name, TString title, Int_t nBins, Double_t lowEdge, Double_t upEdge, TString xtitle);
  void CalculateRatioPlot(TH1F * numer, TH1F * denom, TH1F *& ratioPlot, UInt_t type);
  void DrawSaveTH1Plot(TH1F * histo, TString outputName);
  void PrintTotals();
  void MoveInputToOutDir();

private:
  TString fOutType;
  TString fOutName;
  TFile * fOutRoot;
  TString fInName;
  TFile * fInRoot;
  Bool_t fSaveTotals;
  TCanvas * fTH1Canv;
};
