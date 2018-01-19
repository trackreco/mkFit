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
  void initializeCMSSWEfficiencyTree();
  void initializeCMSSWFakeRateTree();
  void initializeFitTree();
  
  void alignTracks(TrackVec& evt_tracks, TrackExtraVec& evt_extra, bool alignExtra) override;

  void collectFitInfo(const FitVal& tmpfitval, int tkid, int layer) override;

  void resetValidationMaps() override;
  void resetFitBranches();

  void setTrackExtras(Event& ev) override;

  void makeSimTkToRecoTksMaps(Event& ev) override;
  void mapRefTkToRecoTks(const TrackVec& evt_tracks, TrackExtraVec& evt_extras, TkIDToTkIDVecMap& refTkMap);
  void makeSeedTkToRecoTkMaps(Event& ev) override;
  void mapSeedTkToRecoTk(const TrackVec& evt_tracks, const TrackExtraVec& evt_extras, TkIDToTkIDMap& seedTkMap);
  void makeRecoTkToRecoTkMaps(Event& ev) override;
  void makeRecoTkToRecoTkMap(TkIDToTkIDMap& refToPairMap, const TrackVec& reftracks, const TrackExtraVec& refextras, 
			     const TrackVec& pairtracks, const TrackExtraVec& pairextras);
  void makeCMSSWTkToRecoTksMaps(Event& ev) override;
  void makeSeedTkToCMSSWTkMap(Event& ev) override;

  void storeSeedAndMCID(Event& ev);
  void setupCMSSWMatching(const Event & ev, RedTrackVec & reducedCMSSW, LayIdxIDVecMapMap & cmsswHitIDMap);

  int getLastFoundHit(const int trackMCHitID, const int mcTrackID, const Event& ev);

  void fillEfficiencyTree(const Event& ev) override;
  void fillFakeRateTree(const Event& ev) override;
  void fillConfigTree() override;
  void fillCMSSWEfficiencyTree(const Event& ev) override;
  void fillCMSSWFakeRateTree(const Event& ev) override;
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

  // CMSSW to Reco Maps
  TkIDToTkIDVecMap cmsswToBuildMap_;
  TkIDToTkIDVecMap cmsswToFitMap_;

  // Special map for CMSSW tracks to seed track labels --> NOT used for fake rate!!
  TkIDToTkIDMap seedToCmsswMap_;

  // Special map for geting exact CMSSW track that originate build track from seed track through seedIDs
  TkIDToTkIDMap buildToCmsswMap_;
 
  // Special map for associating candidate to fit tracks in CMSSW only
  TkIDToTkIDMap buildToFitMap_;
  TkIDToTkIDMap fitToBuildMap_;

  // Efficiency Tree 
  TTree* efftree_;  
  int   evtID_eff_=0,mcID_eff_=0;
  int   mcmask_seed_eff_=0,mcmask_build_eff_=0,mcmask_fit_eff_=0;
  int   seedID_seed_eff_=0,seedID_build_eff_=0,seedID_fit_eff_=0;

  // for efficiency and duplicate rate plots
  float x_mc_gen_eff_=0.,y_mc_gen_eff_=0.,z_mc_gen_eff_=0.;
  float pt_mc_gen_eff_=0.,phi_mc_gen_eff_=0.,eta_mc_gen_eff_=0.;
  int   nHits_mc_eff_=0,lastlyr_mc_eff_=0;

  // for getting last hit positions track ended up on
  float xhit_seed_eff_=0.,xhit_build_eff_=0.,xhit_fit_eff_=0.;
  float yhit_seed_eff_=0.,yhit_build_eff_=0.,yhit_fit_eff_=0.;
  float zhit_seed_eff_=0.,zhit_build_eff_=0.,zhit_fit_eff_=0.;

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
  int   lastlyr_seed_eff_=0,lastlyr_build_eff_=0,lastlyr_fit_eff_=0;

  // chi2 of tracks
  float hitchi2_seed_eff_=0.,hitchi2_build_eff_=0.,hitchi2_fit_eff_=0.;
  float helixchi2_seed_eff_=0.,helixchi2_build_eff_=0.,helixchi2_fit_eff_=0.;

  // for duplicate track matches
  int   duplmask_seed_eff_=0,duplmask_build_eff_=0,duplmask_fit_eff_=0;
  int   nTkMatches_seed_eff_=0,nTkMatches_build_eff_=0,nTkMatches_fit_eff_=0;

  // Fake Rate tree and variables
  TTree* frtree_;
  int   evtID_FR_=0,seedID_FR_=0;

  int   seedmask_seed_FR_=0,seedmask_build_FR_=0,seedmask_fit_FR_=0;

  // for getting last hit positions track ended up on
  float xhit_seed_FR_=0.,xhit_build_FR_=0.,xhit_fit_FR_=0.;
  float yhit_seed_FR_=0.,yhit_build_FR_=0.,yhit_fit_FR_=0.;
  float zhit_seed_FR_=0.,zhit_build_FR_=0.,zhit_fit_FR_=0.;

  // track state info
  float pt_mc_seed_FR_=0.,pt_mc_build_FR_=0.,pt_mc_fit_FR_=0.;
  float pt_seed_FR_=0.,pt_build_FR_=0.,pt_fit_FR_=0.,ept_seed_FR_=0.,ept_build_FR_=0.,ept_fit_FR_=0.;
  float phi_mc_seed_FR_=0.,phi_mc_build_FR_=0.,phi_mc_fit_FR_=0.;
  float phi_seed_FR_=0.,phi_build_FR_=0.,phi_fit_FR_=0.,ephi_seed_FR_=0.,ephi_build_FR_=0.,ephi_fit_FR_=0.;
  float eta_mc_seed_FR_=0.,eta_mc_build_FR_=0.,eta_mc_fit_FR_=0.;
  float eta_seed_FR_=0.,eta_build_FR_=0.,eta_fit_FR_=0.,eeta_seed_FR_=0.,eeta_build_FR_=0.,eeta_fit_FR_=0.;
    
  int   nHits_seed_FR_=0,nHits_build_FR_=0,nHits_fit_FR_=0;
  int   nHitsMatched_seed_FR_=0,nHitsMatched_build_FR_=0,nHitsMatched_fit_FR_=0;
  float fracHitsMatched_seed_FR_=0,fracHitsMatched_build_FR_=0,fracHitsMatched_fit_FR_=0;
  int   lastlyr_seed_FR_=0,lastlyr_build_FR_=0,lastlyr_fit_FR_=0;

  float hitchi2_seed_FR_=0.,hitchi2_build_FR_=0.,hitchi2_fit_FR_=0.;
 
  int   mcID_seed_FR_=0,mcID_build_FR_=0,mcID_fit_FR_=0;
  int   mcmask_seed_FR_=0,mcmask_build_FR_=0,mcmask_fit_FR_=0;
  int   nHits_mc_seed_FR_=0,nHits_mc_build_FR_=0,nHits_mc_fit_FR_=0;
  int   lastlyr_mc_seed_FR_=0,lastlyr_mc_build_FR_=0,lastlyr_mc_fit_FR_=0;

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

  // CMSSW Efficiency tree
  TTree* cmsswefftree_;
  int   evtID_ceff_=0,cmsswID_ceff_=0,seedID_cmssw_ceff_=0;

  float x_cmssw_ceff_=0.,y_cmssw_ceff_=0.,z_cmssw_ceff_=0.;
  float pt_cmssw_ceff_=0.,phi_cmssw_ceff_=0.,eta_cmssw_ceff_=0.;
  int   nHits_cmssw_ceff_=0,nLayers_cmssw_ceff_=0,lastlyr_cmssw_ceff_=0;

  // build
  int   seedID_build_ceff_=0,mcTrackID_build_ceff_=0;
  int   cmsswmask_build_ceff_=0;

  float pt_build_ceff_=0.,ept_build_ceff_=0.;
  float phi_build_ceff_=0.,ephi_build_ceff_=0.;
  float eta_build_ceff_=0.,eeta_build_ceff_=0.;

  int   nHits_build_ceff_=0,nLayers_build_ceff_=0,nHitsMatched_build_ceff_=0,lastlyr_build_ceff_=0;
  float fracHitsMatched_build_ceff_=0;

  float xhit_build_ceff_=0.,yhit_build_ceff_=0.,zhit_build_ceff_=0.;

  // chi2 of tracks + phi swim
  float hitchi2_build_ceff_=0.,helixchi2_build_ceff_=0.;
  float dphi_build_ceff_=0.;

  int   duplmask_build_ceff_=0,nTkMatches_build_ceff_=0;

  // fit
  int   seedID_fit_ceff_=0,mcTrackID_fit_ceff_=0;
  int   cmsswmask_fit_ceff_=0;

  float pt_fit_ceff_=0.,ept_fit_ceff_=0.;
  float phi_fit_ceff_=0.,ephi_fit_ceff_=0.;
  float eta_fit_ceff_=0.,eeta_fit_ceff_=0.;

  int   nHits_fit_ceff_=0,nLayers_fit_ceff_=0,nHitsMatched_fit_ceff_=0,lastlyr_fit_ceff_=0;
  float fracHitsMatched_fit_ceff_=0;

  float xhit_fit_ceff_=0.,yhit_fit_ceff_=0.,zhit_fit_ceff_=0.;

  // chi2 of tracks + phi swim
  float hitchi2_fit_ceff_=0.,helixchi2_fit_ceff_=0.;
  float dphi_fit_ceff_=0.;

  int   duplmask_fit_ceff_=0,nTkMatches_fit_ceff_=0;

  // CMSSW FakeRate tree
  TTree* cmsswfrtree_;
  int   evtID_cFR_=0,seedID_cFR_=0,mcTrackID_cFR_=0;

  // build info
  int   cmsswID_build_cFR_=0,cmsswmask_build_cFR_=0;

  float pt_build_cFR_=0.,ept_build_cFR_=0.;
  float phi_build_cFR_=0.,ephi_build_cFR_=0.;
  float eta_build_cFR_=0.,eeta_build_cFR_=0.;

  int   nHits_build_cFR_=0,nLayers_build_cFR_=0,nHitsMatched_build_cFR_=0,lastlyr_build_cFR_=0;
  float fracHitsMatched_build_cFR_=0;

  float xhit_build_cFR_=0.,yhit_build_cFR_=0.,zhit_build_cFR_=0.;

  // chi2 of tracks
  float hitchi2_build_cFR_=0.,helixchi2_build_cFR_=0.;
  float dphi_build_cFR_=0.;

  // for duplicate track matches
  int   duplmask_build_cFR_=0,iTkMatches_build_cFR_=0;

  // cmssw info
  int   seedID_cmssw_build_cFR_=0;
  float x_cmssw_build_cFR_=0.,y_cmssw_build_cFR_=0.,z_cmssw_build_cFR_=0.;
  float pt_cmssw_build_cFR_=0.,phi_cmssw_build_cFR_=0.,eta_cmssw_build_cFR_=0.;
  int   nHits_cmssw_build_cFR_=0,nLayers_cmssw_build_cFR_=0,lastlyr_cmssw_build_cFR_=0;

  // fit info
  int   cmsswID_fit_cFR_=0,cmsswmask_fit_cFR_=0;

  float pt_fit_cFR_=0.,ept_fit_cFR_=0.;
  float phi_fit_cFR_=0.,ephi_fit_cFR_=0.;
  float eta_fit_cFR_=0.,eeta_fit_cFR_=0.;

  int   nHits_fit_cFR_=0,nLayers_fit_cFR_=0,nHitsMatched_fit_cFR_=0,lastlyr_fit_cFR_=0;
  float fracHitsMatched_fit_cFR_=0;

  float xhit_fit_cFR_=0.,yhit_fit_cFR_=0.,zhit_fit_cFR_=0.;

  // chi2 of tracks
  float hitchi2_fit_cFR_=0.,helixchi2_fit_cFR_=0.;
  float dphi_fit_cFR_=0.;

  // for duplicate track matches
  int   duplmask_fit_cFR_=0,iTkMatches_fit_cFR_=0;

  // cmssw info
  int   seedID_cmssw_fit_cFR_=0;
  float x_cmssw_fit_cFR_=0.,y_cmssw_fit_cFR_=0.,z_cmssw_fit_cFR_=0.;
  float pt_cmssw_fit_cFR_=0.,phi_cmssw_fit_cFR_=0.,eta_cmssw_fit_cFR_=0.;
  int   nHits_cmssw_fit_cFR_=0,nLayers_cmssw_fit_cFR_=0,lastlyr_cmssw_fit_cFR_=0;

  // Fit tree (for fine tuning z-phi windows and such --> MPlex Only
  TTree* fittree_;
  int   ntotallayers_fit_=0,tkid_fit_=0,evtid_fit_=0;

  static const int   nfvs_ = 24;
  std::vector<float> fvs_[nfvs_];

  // perl -pe 'BEGIN{$i=0} while (s/\(Config::nTotalLayers\)/ = fvs_[$i]/){++$i}' xxin
  std::vector<float> &z_prop_fit_ = fvs_[0],&ez_prop_fit_ = fvs_[1];
  std::vector<float> &z_hit_fit_ = fvs_[2],&ez_hit_fit_ = fvs_[3],&z_sim_fit_ = fvs_[4],&ez_sim_fit_ = fvs_[5];
  std::vector<float> &pphi_prop_fit_ = fvs_[6],&epphi_prop_fit_ = fvs_[7];
  std::vector<float> &pphi_hit_fit_ = fvs_[8],&epphi_hit_fit_ = fvs_[9],&pphi_sim_fit_ = fvs_[10],&epphi_sim_fit_ = fvs_[11];
  std::vector<float> &pt_up_fit_ = fvs_[12],&ept_up_fit_ = fvs_[13],&pt_sim_fit_ = fvs_[14],&ept_sim_fit_ = fvs_[15];
  std::vector<float> &mphi_up_fit_ = fvs_[16],&emphi_up_fit_ = fvs_[17],&mphi_sim_fit_ = fvs_[18],&emphi_sim_fit_ = fvs_[19];
  std::vector<float> &meta_up_fit_ = fvs_[20],&emeta_up_fit_ = fvs_[21],&meta_sim_fit_ = fvs_[22],&emeta_sim_fit_ = fvs_[23];

  std::mutex glock_;
};
#endif
#endif
