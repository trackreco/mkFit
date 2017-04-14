#ifndef _ttreevalidation_
#define _ttreevalidation_

#include "Validation.h"

#ifdef NO_ROOT
class TTreeValidation : public Validation {
public:
  TTreeValidation(std::string) {}
};
#else

#include <unordered_map>
#include <mutex>
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"

// FitVal defined in Validation.h
typedef std::map<int, FitVal> FitValLayMap;
typedef std::unordered_map<int, FitValLayMap> TkIDtoFitValLayMapMap;

class TTreeValidation : public Validation {
public:
  TTreeValidation(std::string fileName);
  ~TTreeValidation();
  
  void initializeEfficiencyTree();
  void initializeFakeRateTree();  
  void initializeConfigTree();
  void initializeFitTree();
  
  void alignTrackExtra(TrackVec& evt_tracks, TrackExtraVec& evt_extra) override;

  void collectFitInfo(const FitVal& tmpfitval, int tkid, int layer) override;

  void resetValidationMaps() override;
  void resetFitBranches();

  void setTrackExtras(Event& ev) override;
  void setTrackCollectionExtras(const TrackVec& evt_tracks, TrackExtraVec& evt_extras, 
				const std::vector<HitVec>& layerHits, const MCHitInfoVec& mcHitInfo);

  void makeSimTkToRecoTksMaps(Event& ev) override;
  void mapSimTkToRecoTks(const TrackVec& evt_tracks, TrackExtraVec& evt_extras, TkIDToTkIDVecMap& simTkMap);
  void makeSeedTkToRecoTkMaps(Event& ev) override;
  void mapSeedTkToRecoTk(const TrackVec& evt_tracks, const TrackExtraVec& evt_extras, TkIDToTkIDMap& seedTkMap);

  int getLastGoodHit(const int trackMCHitID, const int mcTrackID, const Event& ev);

  void fillEfficiencyTree(const Event& ev) override;
  void fillFakeRateTree(const Event& ev) override;
  void fillConfigTree() override;
  void fillFitTree(const Event& ev) override;

  void saveTTrees() override;

 private:
  TFile* f_; // output file!
  
  TkIDtoFitValLayMapMap fitValTkMapMap_; // map used for fit validation in mplex

  // Sim to Reco Maps
  TkIDToTkIDVecMap simToSeedMap_;
  TkIDToTkIDVecMap simToBuildMap_;
  TkIDToTkIDVecMap simToFitMap_;

  // Reco to Reco Maps
  TkIDToTkIDMap seedToBuildMap_;
  TkIDToTkIDMap seedToFitMap_;

  // Efficiency Tree 
  TTree* efftree_;  
  int   evtID_eff_=0,mcID_eff_=0;
  int   mcmask_seed_eff_=0,mcmask_build_eff_=0,mcmask_fit_eff_=0;
  int   seedID_seed_eff_=0,seedID_build_eff_=0,seedID_fit_eff_=0;

  // for efficiency and duplicate rate plots
  float pt_mc_gen_eff_=0.,phi_mc_gen_eff_=0.,eta_mc_gen_eff_=0.;
  int   nHits_mc_eff_=0;

  // for track resolutions / pulls
  float pt_mc_seed_eff_=0.,pt_mc_build_eff_=0.,pt_mc_fit_eff_=0.;
  float pt_seed_eff_=0.,pt_build_eff_=0.,pt_fit_eff_=0.,ept_seed_eff_=0.,ept_build_eff_=0.,ept_fit_eff_=0.;
  float phi_mc_seed_eff_=0.,phi_mc_build_eff_=0.,phi_mc_fit_eff_=0.;
  float phi_seed_eff_=0.,phi_build_eff_=0.,phi_fit_eff_=0.,ephi_seed_eff_=0.,ephi_build_eff_=0.,ephi_fit_eff_=0.;
  float eta_mc_seed_eff_=0.,eta_mc_build_eff_=0.,eta_mc_fit_eff_=0.;
  float eta_seed_eff_=0.,eta_build_eff_=0.,eta_fit_eff_=0.,eeta_seed_eff_=0.,eeta_build_eff_=0.,eeta_fit_eff_=0.;
  
  // for hit countings
  int   nHits_seed_eff_=0,nHits_build_eff_=0,nHits_fit_eff_=0;
  int   nHitsMatched_seed_eff_=0,nHitsMatched_build_eff_=0,nHitsMatched_fit_eff_=0;
  float fracHitsMatched_seed_eff_=0,fracHitsMatched_build_eff_=0,fracHitsMatched_fit_eff_=0;

  // chi2 of tracks
  float hitchi2_seed_eff_=0.,hitchi2_build_eff_=0.,hitchi2_fit_eff_=0.;
  float helixchi2_seed_eff_=0.,helixchi2_build_eff_=0.,helixchi2_fit_eff_=0.;

  // for duplicate track matches
  int   duplmask_seed_eff_=0,duplmask_build_eff_=0,duplmask_fit_eff_=0;
  int   nTkMatches_seed_eff_=0,nTkMatches_build_eff_=0,nTkMatches_fit_eff_=0;

  // Fake Rate tree and variables
  TTree* fakeratetree_;
  int   evtID_FR_=0,seedID_FR_=0;

  int   seedmask_seed_FR_=0,seedmask_build_FR_=0,seedmask_fit_FR_=0;
  float pt_mc_seed_FR_=0.,pt_mc_build_FR_=0.,pt_mc_fit_FR_=0.;
  float pt_seed_FR_=0.,pt_build_FR_=0.,pt_fit_FR_=0.,ept_seed_FR_=0.,ept_build_FR_=0.,ept_fit_FR_=0.;
  float phi_mc_seed_FR_=0.,phi_mc_build_FR_=0.,phi_mc_fit_FR_=0.;
  float phi_seed_FR_=0.,phi_build_FR_=0.,phi_fit_FR_=0.,ephi_seed_FR_=0.,ephi_build_FR_=0.,ephi_fit_FR_=0.;
  float eta_mc_seed_FR_=0.,eta_mc_build_FR_=0.,eta_mc_fit_FR_=0.;
  float eta_seed_FR_=0.,eta_build_FR_=0.,eta_fit_FR_=0.,eeta_seed_FR_=0.,eeta_build_FR_=0.,eeta_fit_FR_=0.;
    
  int   nHits_seed_FR_=0,nHits_build_FR_=0,nHits_fit_FR_=0;
  int   nHitsMatched_seed_FR_=0,nHitsMatched_build_FR_=0,nHitsMatched_fit_FR_=0;
  float fracHitsMatched_seed_FR_=0,fracHitsMatched_build_FR_=0,fracHitsMatched_fit_FR_=0;

  float hitchi2_seed_FR_=0.,hitchi2_build_FR_=0.,hitchi2_fit_FR_=0.;
 
  int   mcID_seed_FR_=0,mcID_build_FR_=0,mcID_fit_FR_=0;
  int   mcmask_seed_FR_=0,mcmask_build_FR_=0,mcmask_fit_FR_=0;
  int   nHits_mc_seed_FR_=0,nHits_mc_build_FR_=0,nHits_mc_fit_FR_=0;

  float helixchi2_seed_FR_=0.,helixchi2_build_FR_=0.,helixchi2_fit_FR_=0.;

  int   duplmask_seed_FR_=0,duplmask_build_FR_=0,duplmask_fit_FR_=0;
  int   iTkMatches_seed_FR_=0,iTkMatches_build_FR_=0,iTkMatches_fit_FR_=0;

  // Configuration tree
  TTree* configtree_;
  int   Ntracks_=0,Nevents_=0;
  int   nLayers_=0;
  float fRadialSpacing_=0.,fRadialExtent_=0.,fInnerSensorSize_=0.,fOuterSensorSize_=0.;
  float fEtaDet_=0.,fPhiFactor_=0.;
  int   nPhiPart_=0,nEtaPart_=0;
  int   nlayers_per_seed_=0,maxCand_=0;
  float chi2Cut_=0.,nSigma_=0.,minDPhi_=0.,maxDPhi_=0.,minDEta_=0.,maxDEta_=0.;
  float beamspotX_=0.,beamspotY_=0.,beamspotZ_=0.;
  float minSimPt_=0.,maxSimPt_=0.;
  float hitposerrXY_=0.,hitposerrZ_=0.,hitposerrR_=0.;
  float varXY_=0.,varZ_=0.;
  int   nTotHit_=0;
  float ptinverr049_=0.,phierr049_=0.,thetaerr049_=0.,ptinverr012_=0.,phierr012_=0.,thetaerr012_=0.;

  // Fit tree (for fine tuning z-phi windows and such --> MPlex Only
  TTree* fittree_;
  int   ntotallayers_fit_=0,tkid_fit_=0,evtid_fit_=0;
  float z_prop_fit_[Config::nTotalLayers],ez_prop_fit_[Config::nTotalLayers];
  float z_hit_fit_[Config::nTotalLayers],ez_hit_fit_[Config::nTotalLayers],z_sim_fit_[Config::nTotalLayers],ez_sim_fit_[Config::nTotalLayers];
  float pphi_prop_fit_[Config::nTotalLayers],epphi_prop_fit_[Config::nTotalLayers];
  float pphi_hit_fit_[Config::nTotalLayers],epphi_hit_fit_[Config::nTotalLayers],pphi_sim_fit_[Config::nTotalLayers],epphi_sim_fit_[Config::nTotalLayers];
  float pt_up_fit_[Config::nTotalLayers],ept_up_fit_[Config::nTotalLayers],pt_sim_fit_[Config::nTotalLayers],ept_sim_fit_[Config::nTotalLayers];
  float mphi_up_fit_[Config::nTotalLayers],emphi_up_fit_[Config::nTotalLayers],mphi_sim_fit_[Config::nTotalLayers],emphi_sim_fit_[Config::nTotalLayers];
  float meta_up_fit_[Config::nTotalLayers],emeta_up_fit_[Config::nTotalLayers],meta_sim_fit_[Config::nTotalLayers],emeta_sim_fit_[Config::nTotalLayers];

  std::mutex glock_;
};
#endif
#endif
