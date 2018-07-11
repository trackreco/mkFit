#ifndef _PlotValidation_
#define _PlotValidation_

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TString.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

typedef std::vector<Float_t>  FltVec;
typedef std::vector<FltVec>   FltVecVec;
typedef std::vector<Double_t> DblVec;
typedef std::vector<DblVec>   DblVecVec;
typedef std::vector<Int_t>    IntVec;
typedef std::vector<TString>  TStrVec;

typedef std::vector<TBranch*>    TBrRefVec;
typedef std::vector<TBrRefVec>   TBrRefVecVec;
typedef std::vector<TDirectory*> TDirRefVec;

typedef std::map<TString,TH1F*>        TH1FRefMap;
typedef std::map<TString,TEfficiency*> TEffRefMap;

struct EffStruct
{
  EffStruct(){}
  ~EffStruct(){}

  Float_t passed_;
  Float_t total_;

  Float_t eff_;
  Float_t elow_;
  Float_t eup_;
};

class PlotValidation
{
public:
  PlotValidation(const TString & inName, const TString & outName, const Bool_t cmsswComp,
		 const Bool_t mvInput, const Bool_t saveAs, const TString & outType);
  ~PlotValidation();
  
  // setup functions
  void SetupStyle();
  void SetupBins();
  void SetupVariableBins(const std::string & s_bins, DblVec & bins);
  void SetupFixedBins(const Int_t nBins, const Double_t low, const Double_t high, DblVec & bins);

  // main call
  void Validation();
  void PlotEffTree();
  void PlotFRTree();
  void PrintTotals();

  // output functions
  template <typename T>
  void DrawWriteSavePlot(T *& plot, TDirectory *& subdir, const TString & subdirname, const TString & option);

  // helper functions
  void MakeOutDir(const TString & outdirname);
  void GetTotalEfficiency(const TEfficiency * eff, EffStruct & effs);
  TDirectory * MakeSubDirs(const TString & subdirname);
  void MoveInput();

private:
  // input+output config
  const TString fInName;
  const Bool_t  fCmsswComp;
  const Bool_t  fMvInput;
  const Bool_t  fSaveAs;
  const TString fOutType;

  // main input 
  TFile * fInRoot;

  // binning for rate plots
  DblVec fPtBins;
  DblVec fEtaBins;
  DblVec fPhiBins;

  // binning for hit hists
  DblVec fNHitsBins;
  DblVec fFracHitsBins;

  // binning for diff hists
  DblVec fDNHitsBins;
  DblVec fDInvPtBins;
  DblVec fDPhiBins;
  DblVec fDEtaBins;

  // output variables
  TString fOutName;
  TFile * fOutRoot;
};

#endif
