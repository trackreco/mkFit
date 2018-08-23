#ifndef _StackValidation_
#define _StackValidation_

#include "Common.hh"

#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

struct RateOpts
{
  RateOpts() {}
  RateOpts(const TString & dir, const TString & sORr, const TString & rate)
    : dir(dir), sORr(sORr), rate(rate) {}

  TString dir;
  TString sORr; // sim or reco 
  TString rate;
};
typedef std::vector<RateOpts> ROVec;

namespace
{
  ROVec rates;
  void setupRates(const Bool_t cmsswComp)
  {
    const TString ref = (cmsswComp?"cmssw":"sim");

    rates.emplace_back("efficiency",ref,"eff");
    rates.emplace_back("inefficiency",ref,"ineff_brl");
    rates.emplace_back("inefficiency",ref,"ineff_trans");
    rates.emplace_back("inefficiency",ref,"ineff_ec");
    rates.emplace_back("fakerate","reco","fr");
    rates.emplace_back("duplicaterate",ref,"dr");
    
    if (cmsswComp) 
    {
      for (UInt_t i = 0; i < rates.size(); i++) rates[i].dir += "_cmssw";
    }
  }

  std::vector<Float_t> ptcuts;
  void setupPtCuts()
  {
    ptcuts.emplace_back(0.f);
    ptcuts.emplace_back(0.9f);
    ptcuts.emplace_back(2.f);
  }
};

class StackValidation
{
public:
  StackValidation(const TString & label, const TString & extra, const Bool_t cmsswComp, const TString & suite);
  ~StackValidation();
  void MakeValidationStacks();
  void MakeRatioStacks(const TString & trk);
  void MakeCMSSWKinematicDiffStacks(const TString & trk);

private:
  const TString label;
  const TString extra;
  const Bool_t cmsswComp;
  const TString suite;

  // legend height
  Double_t y1;
  Double_t y2; 

  std::vector<TFile*> files;
};

#endif
