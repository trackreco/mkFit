#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

#include <vector>
#include <iostream>
#include <cmath>

typedef std::vector<TString> TStrVec;
typedef std::vector<UInt_t>  UIntVec;
typedef std::vector<Int_t>   IntVec;
typedef std::vector<Float_t> FltVec;
typedef std::vector<FltVec>  FltVecVec;
typedef std::vector<std::vector<TH1F *> > TH1FRefVecVec;

class PlotValidation
{
public:
  PlotValidation(TString inName, TString outName, TString outType);
  ~PlotValidation();
  void Validation(Bool_t mvInput = false);
  void PlotResPull();
  void PlotEfficiency();
  void PlotFakeRate();
  void ComputeResPull(const Float_t& mc_val, const Float_t& var_val, const Float_t& var_err, Float_t var_out[]);
  void ComputeRatioPlot(const TH1F * numer, const TH1F * denom, TH1F *& ratioPlot);
  void DrawSaveTH1Plot(TH1F * histo, const TString plotName);
  void DrawFitSaveTH1Plot(TH1F * histo, const TString plotName, const Float_t fitrange);
  void MoveInput();

private:
  TString fInName;
  TFile * fInRoot;
  TString fOutType;
  TString fOutName;
  TFile * fOutRoot;
  TCanvas * fTH1Canv;
};
