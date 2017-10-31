#ifndef _PlotValidation_
#define _PlotValidation_

#include "TFile.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TColor.h"
#include "TLegend.h"

#include <vector>
#include <iostream>
#include <iomanip>
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
typedef std::vector<TH1FRefVecVecVec> TH1FRefVecVecVecVec;

typedef std::vector<TEfficiency *> TEffRefVec;
typedef std::vector<TEffRefVec>    TEffRefVecVec;
typedef std::vector<TEffRefVecVec> TEffRefVecVecVec;
typedef std::vector<TEffRefVecVecVec> TEffRefVecVecVecVec;

typedef std::vector<TString> TStrVec;

struct EffStruct
{
  Float_t passed_;
  Float_t total_;

  Float_t eff_;
  Float_t elow_;
  Float_t eup_;
};

class PlotValidation
{
public:
  PlotValidation(TString inName, TString outName, Bool_t computePulls, Bool_t cmsswComp,
		 Bool_t mvInput, Bool_t saveAs, TString outType);
  ~PlotValidation();
  void Validation();

  void PlotEfficiency();
  void PlotInefficiencyVsGeom();
  void PlotFakeRate();
  void PlotDuplicateRate();
  void PlotNHits();
  void PlotMomResolutionPull();
  void PlotCMSSWKinematicDiffs(); 

  void PrintTotals();

  void MakeSubDirectory(const TString subdirname);

  void ComputeResidual      (const Float_t mcvar_val, const Float_t recovar_val, Float_t & var_out);
  void ComputeResolutionPull(const Float_t mcvar_val, const Float_t recovar_val, const Float_t recovar_err, FltVec & var_out);

  void GetTotalEfficiency(const TEfficiency * eff, EffStruct & effs);
    
  void DrawWriteSaveTEffPlot   (TDirectory *& subdir, TEfficiency *& eff, const TString subdirname, const TString plotName);
  void DrawWriteSaveTH1FPlot   (TDirectory *& subdir, TH1F *& histo, const TString subdirname, const TString plotName);
  void DrawWriteFitSaveTH1FPlot(TDirectory *& subdir, TH1F *& histo, const TString subdirname, const TString plotName, const Float_t fitrange);

  void MoveInput();

private:
  const TString fInName;
        TFile * fInRoot;
  const Bool_t  fComputePulls;
  const Bool_t  fCmsswComp;
  const Bool_t  fMvInput;
  const Bool_t  fSaveAs;
  const TString fOutType;
  TString fOutName;
  TFile * fOutRoot;
  TCanvas * fTEffCanv;
  TCanvas * fTH1Canv;

  // color base
  std::vector<Color_t> fColors;
  UInt_t fColorSize;
};

#endif
