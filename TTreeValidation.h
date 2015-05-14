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

typedef std::map<unsigned int,TrackVecRef> simToTkMap;

class TTreeValidation : public Validation {
public:
  TTreeValidation(std::string fileName);

  void fillBuildTree(const unsigned int layer, const unsigned int branches, const unsigned int cands) override;
  void makeSimToTkMaps(TrackVec& evt_seed_tracks, TrackVec& evt_build_tracks, TrackVec& evt_fit_tracks) override;
  void mapSimToTks(TrackVec& evt_tracks, simToTkMap& simTkMap_);
  void fillEffTree(const TrackVec& evt_sim_tracks, const unsigned int & ievt) override;
  void fillFakeTrees(const unsigned int & ievt) override;
  void fillFakeSeedTree(const unsigned int & ievt);
  void fillFakeBuildTree(const unsigned int & ievt);
  void fillFakeFitTree(const unsigned int & ievt);
  void fillConfigTree(const unsigned int & ntracks, const unsigned int & nevts, const std::vector<double> & ticks);
  void saveTTrees() override;

  simToTkMap simToSeedMap_;
  simToTkMap simToBuildMap_;
  simToTkMap simToFitMap_;

  TFile* f_;

  TTree* configtree_;
  float        simtime_=0.,segtime_=0.,seedtime_=0.,buildtime_=0.,fittime_=0.,hlvtime_=0.;
  unsigned int Ntracks_=0,Nevents_=0;
  unsigned int nPhiPart_=0,nPhiFactor_=0,nEtaPart_=0;
  float        etaDet_=0.;
  unsigned int nlayers_per_seed_=0,maxCand_=0;
  float        chi2Cut_=0.,nSigma_=0.,minDPhi_=0.,maxDPhi_=0.,minDEta_=0.,maxDEta_=0.;

  TTree* tree_br_;
  unsigned int layer_=0,branches_=0,cands_=0;

  TTree* efftree_;  
  //  std::vector<float> pos_mc_;  
  float pt_mc_eff_=0.,pt_seed_eff_=0.,pt_build_eff_=0.,pt_fit_eff_=0.,ept_seed_eff_=0.,ept_build_eff_=0.,ept_fit_eff_=0.;
  float pz_mc_eff_=0.,pz_seed_eff_=0.,pz_build_eff_=0.,pz_fit_eff_=0.,epz_seed_eff_=0.,epz_build_eff_=0.,epz_fit_eff_=0.;
  float phi_mc_eff_=0.,phi_seed_eff_=0.,phi_build_eff_=0.,phi_fit_eff_=0.,ephi_seed_eff_=0.,ephi_build_eff_=0.,ephi_fit_eff_=0.;
  float eta_mc_eff_=0.,eta_seed_eff_=0.,eta_build_eff_=0.,eta_fit_eff_=0.,eeta_seed_eff_=0.,eeta_build_eff_=0.,eeta_fit_eff_=0.;
  int   nHits_mc_eff_=0,nHits_seed_eff_=0,nHits_build_eff_=0,nHits_fit_eff_=0;
  float chi2_seed_eff_=0.,chi2_build_eff_=0.,chi2_fit_eff_=0.;
  int   nHitsMatched_seed_eff_=0,nHitsMatched_build_eff_=0,nHitsMatched_fit_eff_=0;
  int   nDup_seed_eff_=0,nDup_build_eff_=0,nDup_fit_eff_=0;
  int   mask_seed_eff_=0,mask_build_eff_=0,mask_fit_eff_=0;
  int   evt_mc_eff_=0,evt_seed_eff_=0,evt_build_eff_=0,evt_fit_eff_=0;

  TTree* fakeseedtree_;
  TTree* fakebuildtree_;
  TTree* fakefittree_;  
  bool  mask_seed_real_=0,mask_build_real_=0,mask_fit_real_=0;
  int   mask_seed_duplicate_=0,mask_build_duplicate_=0,mask_fit_duplicate_=0;
  float pt_seed_fake_=0.,pt_build_fake_=0.,pt_fit_fake_=0.;
  float pz_seed_fake_=0.,pz_build_fake_=0.,pz_fit_fake_=0.;
  float phi_seed_fake_=0.,phi_build_fake_=0.,phi_fit_fake_=0.;
  float eta_seed_fake_=0.,eta_build_fake_=0.,eta_fit_fake_=0.;
  int   nHits_seed_fake_=0,nHits_build_fake_=0,nHits_fit_fake_=0;
  int   nHitsMatched_seed_fake_=0,nHitsMatched_build_fake_=0,nHitsMatched_fit_fake_=0;
  float chi2_seed_fake_=0.,chi2_build_fake_=0.,chi2_fit_fake_=0.;
  int   evt_seed_fake_=0,evt_build_fake_=0,evt_fit_fake_=0;

  std::mutex glock_;
};
#endif
#endif
