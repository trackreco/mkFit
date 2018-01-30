#ifndef _StackValidation_
#define _StackValidation_

#include "plotting/Common.hh"

#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

struct RateOpts
{
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

    rates.push_back({"efficiency",ref,"eff"});
    rates.push_back({"inefficiency",ref,"ineff_barrel"});
    rates.push_back({"inefficiency",ref,"ineff_endcap"});
    rates.push_back({"fakerate","reco","fr"});
    rates.push_back({"duplicaterate",ref,"dr"});
    
    if (cmsswComp) 
    {
      for (UInt_t i = 0; i < rates.size(); i++) rates[i].dir += "_cmssw";
    }
  }

  std::vector<Float_t> ptcuts;
  void setupPtCuts()
  {
    ptcuts.push_back(0.f);
    ptcuts.push_back(0.9f);
    ptcuts.push_back(2.f);
  }
};

class StackValidation
{
public:
  StackValidation(const TString & label, const TString & extra, const Bool_t cmsswComp);
  ~StackValidation();
  void MakeValidationStacks();
  void MakeRatioStacks(const TString & trk);
  void MakeCMSSWKinematicDiffStacks(const TString & trk);

private:
  const TString label;
  const TString extra;
  const Bool_t cmsswComp;

  std::vector<TFile*> files;
};

#endif
