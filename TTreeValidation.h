#ifndef _ttreevalidation_
#define _ttreevalidation_

#include "Validation.h"

#ifdef NO_ROOT
class TTreeValidation : public Validation {
public:
  TTreeValidation(std::string) {}
};
#else

#include <map>
#include <mutex>
#include "TFile.h"
#include "TTree.h"

typedef std::map<unsigned int,TrackVec> simToTkMap;

class TTreeValidation : public Validation {
public:
  TTreeValidation(std::string fileName);

  void fillBuildTree(const unsigned int layer, const unsigned int branches, const unsigned int cands) override;
  void makeSimToTkMaps(TrackVec& evt_seed_tracks, TrackVec& evt_build_tracks, TrackVec& evt_fit_tracks) override;
  void mapSimToTks(TrackVec& evt_tracks, simToTkMap& simTkMap_);
  void fillEffTree(const TrackVec& evt_sim_tracks, const TrackVec& evt_seed_tracks, const TrackVec& evt_build_tracks, const TrackVec& evt_fit_tracks, const unsigned int & ievt) override;
  void saveTTrees() override;

  simToTkMap simToSeedMap_;
  simToTkMap simToBuildMap_;
  simToTkMap simToFitMap_;

  TFile* f_;
  TTree* tree_br_;
  unsigned int layer_=0,branches_=0,cands_=0;

  TTree* efftree_;  
  float pt_mc_=0.,pt_seed_=0.,pt_build_=0.,pt_fit_=0.,ept_seed_=0.,ept_build_=0.,ept_fit_=0.;
  float pz_mc_=0.,pz_seed_=0.,pz_build_=0.,pz_fit_=0.,epz_seed_=0.,epz_build_=0.,epz_fit_=0.;
  float phi_mc_=0.,phi_seed_=0.,phi_build_=0.,phi_fit_=0.,ephi_seed_=0.,ephi_build_=0.,ephi_fit_=0.;
  float eta_mc_=0.,eta_seed_=0.,eta_build_=0.,eta_fit_=0.,eeta_seed_=0.,eeta_build_=0.,eeta_fit_=0.;
  int   nHits_mc_=0,nHits_seed_=0,nHits_build_=0,nHits_fit_=0;
  float chi2_seed_=0.,chi2_build_=0.,chi2_fit_=0.;
  int   nHitsMatched_seed_=0,nHitsMatched_build_=0,nHitsMatched_fit_=0;
  int   evt_mc_=0,evt_seed_=0,evt_build_=0,evt_fit_=0;

  std::mutex glock_;
};
#endif
#endif
