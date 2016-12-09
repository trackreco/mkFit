#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TColor.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

typedef std::vector<Float_t>      FltVec;
typedef std::vector<FltVec>       FltVecVec;
typedef std::vector<FltVecVec>    FltVecVecVec;
typedef std::vector<Float_t> *    FltVecRef;
typedef std::vector<FltVecRef>    FltVecRefVec;
typedef std::vector<FltVecRefVec> FltVecRefVecVec;

typedef std::vector<Int_t>     IntVec;
typedef std::vector<IntVec>    IntVecVec;
typedef std::vector<Int_t> *   IntVecRef;
typedef std::vector<IntVecRef> IntVecRefVec;

typedef std::vector<TH1F *>        TH1FRefVec;
typedef std::vector<TH1FRefVec>    TH1FRefVecVec;
typedef std::vector<TH1FRefVecVec> TH1FRefVecVecVec;

typedef std::vector<TString> TStrVec;

class PlotValidation
{
public:
  PlotValidation(TString inName, TString outName, 
		 Bool_t mvInput, Bool_t fullVal,
		 Bool_t saveAs, TString outType);
  ~PlotValidation();
  void Validation();

  void PlotEfficiency();
  void PlotFakeRate();
  void PlotDuplicateRate();
  void PlotNHits();
  void PlotTiming();
  void PlotMomResolutionPull();

  void PlotSegment();
  void PlotBranching();
  void PlotSimGeo();
  void PlotPosResolutionPull();
  void PlotCFResidual();
  void PlotCFResolutionPull();

  void PrintTotals();

  void MakeSubDirectory(const TString subdirname);

  void ComputeResidual      (const Float_t mcvar_val, const Float_t recovar_val, Float_t & var_out);
  void ComputeResolutionPull(const Float_t mcvar_val, const Float_t recovar_val, const Float_t recovar_err, FltVec & var_out);

  void ComputeRatioPlot(const TH1F * numer, const TH1F * denom, TH1F *& ratioPlot, Bool_t subone = false);

  void ZeroSuppressPlot(TH1F *& histo);

  void WriteTH2FPlot        (TDirectory *& subdir, TH2F *& hist);
  void DrawWriteSaveTH2FPlot(TDirectory *& subdir, TH2F *& histo, const TString subdirname, const TString plotName);

  void WriteTH1FPlot           (TDirectory *& subdir, TH1F *& histo);
  void DrawWriteSaveTH1FPlot   (TDirectory *& subdir, TH1F *& histo, const TString subdirname, const TString plotName, const Bool_t zeroSupLin);
  void DrawWriteFitSaveTH1FPlot(TDirectory *& subdir, TH1F *& histo, const TString subdirname, const TString plotName, const Float_t fitrange);

  void MoveInput();

private:
  TString fInName;
  TFile * fInRoot;
  Bool_t  fMvInput;
  Bool_t  fFullVal;
  Bool_t  fSaveAs;
  TString fOutType;
  TString fOutName;
  TFile * fOutRoot;
  TCanvas * fTH1Canv;
  TCanvas * fTH2Canv;

  // color base
  std::vector<Color_t> fColors;
  UInt_t fColorSize;

};
