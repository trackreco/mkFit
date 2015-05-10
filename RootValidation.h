#ifndef _rootvalidation_
#define _rootvalidation_

#include "Validation.h"

#ifdef NO_ROOT
class RootValidation : public Validation {
public:
  RootValidation(std::string) {}
};
#else

#include <map>
#include <mutex>
#include "TFile.h"
#include "TTree.h"

typedef std::map<unsigned int,TrackVec> simToTkMap;

class RootValidation : public Validation {
public:
  RootValidation(std::string fileName);

  void fillBuildTree(unsigned int layer, unsigned int branches, unsigned int cands) override;
  void fillEffTree(TrackVec& evt_sim_tracks, TrackVec& evt_seed_tracks, TrackVec& evt_build_tracks, TrackVec& evt_fit_tracks) override;
  void saveTTrees() override;
  
  void mapSimToTks(TrackVec& evt_tracks, simToTkMap& simTkMap);

  TFile* f_;
  TTree* tree_br_;
  unsigned int layer_=0,branches_=0,cands_=0;

  TTree* efftree_;  
  float pt_mc_=0.,pt_seed_=0.,pt_build_=0.,pt_fit_=0.,ept_seed_=0.,ept_build_=0.,ept_fit_=0.;
  float pz_mc_=0.,pz_seed_=0.,pz_build_=0.,pz_fit_=0.,epz_seed_=0.,epz_build_=0.,epz_fit_=0.;
  float phi_mc_=0.,phi_seed_=0.,phi_build_=0.,phi_fit_=0.,ephi_seed_=0.,ephi_build_=0.,ephi_fit_=0.;
  float eta_mc_=0.,eta_seed_=0.,eta_build_=0.,eta_fit_=0.,eeta_seed_=0.,eeta_build_=0.,eeta_fit_=0.;
  int   nHits_mc_=0,nHits_seed_=0,nHits_build_=0,nHits_fit_=0;
  float chi2_mc_=0.,chi2_seed_=0.,chi2_build_=0.,chi2_fit_=0.;
  
  std::mutex glock_;
};
#endif
#endif
