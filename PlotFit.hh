#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

#include <iostream>
#include <cmath>

class PlotFit
{
public:
  PlotFit(TString inName, TString outName, TString outType);
  ~PlotFit();
  void PlotAllHistos();
  void PlotPosResPull();
  void PlotOverFlow();
  void PlotPtResPull();
  void PlotGeo();
  void ComputeResPull(const Float_t& init_val, const Float_t& step_val, const Float_t& step_err, Float_t step_out[]);
  void DrawFitSaveTH1Plot(TH1F * histo, Float_t fitRange, TString plotName);
  void DrawSaveTH1Plot(TH1F * histo, TString plotName);
  void DrawSaveTH2Plot(TH2F * histo, TString plotName, Int_t pixelRange[]);

private:
  TString fOutType;
  TString fOutName;
  TFile * fOutRoot;
  TFile * fInRoot;
  TCanvas * fTH1Canv;
  TCanvas * fTH2Canv;
};
