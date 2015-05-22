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
#include "TROOT.h"

typedef std::map<unsigned int,TrackVecRef> simToTksMap;
typedef std::map<unsigned int,const Track*> seedToTkMap;

class TTreeValidation : public Validation {
public:
  TTreeValidation(std::string fileName);

  void fillBuildTree(const unsigned int layer, const unsigned int branches, const unsigned int cands) override;

  void makeSimToTksMaps(TrackVec& evt_seed_tracks, TrackVec& evt_build_tracks, TrackVec& evt_fit_tracks) override;
  void mapSimToTks(TrackVec& evt_tracks, simToTksMap& simTkMap);
  void fillEffTree(const TrackVec& evt_sim_tracks, const unsigned int & ievt) override;

  void makeSeedToTkMaps(const TrackVec& evt_build_tracks, const TrackVec& evt_fit_tracks) override;
  void mapSeedToTk(const TrackVec& evt_tracks, seedToTkMap& seedTkMap);
  void fillFakeRateTree(const TrackVec& evt_sim_tracks, const TrackVec& evt_seed_tracks, const unsigned int & ievt) override;

  void fillChi2Tree(const TrackVec & evt_sim_tracks, const unsigned int& ievt) override;

  void fillConfigTree(const unsigned int & ntracks, const unsigned int & nevts, const std::vector<double> & ticks);

  void saveTTrees() override;

  TFile* f_;

  // build branching tree
  TTree* tree_br_;
  unsigned int layer_=0,branches_=0,cands_=0;

  // efficiency trees and variables
  simToTksMap simToSeedMap_;
  simToTksMap simToBuildMap_;
  simToTksMap simToFitMap_;

  TTree* efftree_;  
  unsigned int evtID_eff_=0,mcID_eff_=0;
  int   mcmask_seed_eff_=0,mcmask_build_eff_=0,mcmask_fit_eff_=0;

  std::vector<float> x_hit_mc_eff_;
  std::vector<float> y_hit_mc_eff_,z_hit_mc_eff_;
  std::vector<float> x_inithit_mc_eff_,y_inithit_mc_eff_,z_inithit_mc_eff_;
  float x_vrx_mc_eff_=0.,y_vrx_mc_eff_=0.,z_vrx_mc_eff_=0.;

  float pt_mc_eff_=0.,pt_seed_eff_=0.,pt_build_eff_=0.,pt_fit_eff_=0.,ept_seed_eff_=0.,ept_build_eff_=0.,ept_fit_eff_=0.;
  float px_mc_eff_=0.,px_seed_eff_=0.,px_build_eff_=0.,px_fit_eff_=0.,epx_seed_eff_=0.,epx_build_eff_=0.,epx_fit_eff_=0.;
  float py_mc_eff_=0.,py_seed_eff_=0.,py_build_eff_=0.,py_fit_eff_=0.,epy_seed_eff_=0.,epy_build_eff_=0.,epy_fit_eff_=0.;
  float pz_mc_eff_=0.,pz_seed_eff_=0.,pz_build_eff_=0.,pz_fit_eff_=0.,epz_seed_eff_=0.,epz_build_eff_=0.,epz_fit_eff_=0.;
  float phi_mc_eff_=0.,phi_seed_eff_=0.,phi_build_eff_=0.,phi_fit_eff_=0.,ephi_seed_eff_=0.,ephi_build_eff_=0.,ephi_fit_eff_=0.;
  float eta_mc_eff_=0.,eta_seed_eff_=0.,eta_build_eff_=0.,eta_fit_eff_=0.,eeta_seed_eff_=0.,eeta_build_eff_=0.,eeta_fit_eff_=0.;

  float chi2_seed_eff_=0.,chi2_build_eff_=0.,chi2_fit_eff_=0.;
  int   nHits_mc_eff_=0,nHits_seed_eff_=0,nHits_build_eff_=0,nHits_fit_eff_=0;
  int   nHitsMatched_seed_eff_=0,nHitsMatched_build_eff_=0,nHitsMatched_fit_eff_=0;
  int   duplmask_seed_eff_=0,duplmask_build_eff_=0,duplmask_fit_eff_=0;
  int   nDup_seed_eff_=0,nDup_build_eff_=0,nDup_fit_eff_=0;

  // fake rate tree and variables
  seedToTkMap seedToBuildMap_;
  seedToTkMap seedToFitMap_;
  
  TTree* fakeratetree_;
  unsigned int evtID_FR_=0,seedID_FR_=0;

  int   seedmask_build_FR_=0,seedmask_fit_FR_=0;
  float pt_seed_FR_=0.,pt_build_FR_=0.,pt_fit_FR_=0.,ept_seed_FR_=0.,ept_build_FR_=0.,ept_fit_FR_=0.;
  float pz_seed_FR_=0.,pz_build_FR_=0.,pz_fit_FR_=0.,epz_seed_FR_=0.,epz_build_FR_=0.,epz_fit_FR_=0.;
  float phi_seed_FR_=0.,phi_build_FR_=0.,phi_fit_FR_=0.,ephi_seed_FR_=0.,ephi_build_FR_=0.,ephi_fit_FR_=0.;
  float eta_seed_FR_=0.,eta_build_FR_=0.,eta_fit_FR_=0.,eeta_seed_FR_=0.,eeta_build_FR_=0.,eeta_fit_FR_=0.;

  float chi2_seed_FR_=0.,chi2_build_FR_=0.,chi2_fit_FR_=0.;
  int   nHits_seed_FR_=0,nHits_build_FR_=0,nHits_fit_FR_=0;

  int   mcID_seed_FR_=0,mcID_build_FR_=0,mcID_fit_FR_=0;
  int   mcmask_seed_FR_=0,mcmask_build_FR_=0,mcmask_fit_FR_=0;
  float pt_mc_seed_FR_=0.,pt_mc_build_FR_=0.,pt_mc_fit_FR_=0.;
  float pz_mc_seed_FR_=0.,pz_mc_build_FR_=0.,pz_mc_fit_FR_=0.;
  float phi_mc_seed_FR_=0.,phi_mc_build_FR_=0.,phi_mc_fit_FR_=0.;
  float eta_mc_seed_FR_=0.,eta_mc_build_FR_=0.,eta_mc_fit_FR_=0.;

  int   nHitsMatched_seed_FR_=0,nHitsMatched_build_FR_=0,nHitsMatched_fit_FR_=0;
  int   nHits_mc_seed_FR_=0,nHits_mc_build_FR_=0,nHits_mc_fit_FR_=0;

  int   duplmask_seed_FR_=0,duplmask_build_FR_=0,duplmask_fit_FR_=0;
  int   iDup_seed_FR_=0,iDup_build_FR_=0,iDup_fit_FR_=0;

  // configuration tree
  TTree* configtree_;
  float        simtime_=0.,segtime_=0.,seedtime_=0.,buildtime_=0.,fittime_=0.,hlvtime_=0.;
  unsigned int Ntracks_=0,Nevents_=0;
  unsigned int nPhiPart_=0,nPhiFactor_=0,nEtaPart_=0;
  float        etaDet_=0.;
  unsigned int nlayers_per_seed_=0,maxCand_=0;
  float        chi2Cut_=0.,nSigma_=0.,minDPhi_=0.,maxDPhi_=0.,minDEta_=0.,maxDEta_=0.;

  // helix chi2 tree
  TTree* helixchi2tree_;
  unsigned int evtID_chi2_=0,mcID_chi2_=0;
  float        helixchi2_seed_=0,helixchi2_build_=0,helixchi2_r0_fit_=0,helixchi2_r1_fit_=0;

  std::mutex glock_;
};
#endif
#endif
