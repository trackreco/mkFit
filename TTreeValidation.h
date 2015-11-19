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

// branching valiation typedefs
struct BranchVal
{
public:
  BranchVal() {}
  float nSigmaDeta;
  float nSigmaDphi;

  unsigned int etaBinMinus;
  unsigned int etaBinPlus;
  unsigned int phiBinMinus;
  unsigned int phiBinPlus;

  std::vector<unsigned int> cand_hit_indices;
  std::vector<unsigned int> branch_hit_indices; // size of was just branches
};

typedef std::vector<BranchVal> BranchValVec;
typedef std::unordered_map<unsigned int, BranchValVec> BranchValVecLayMap;
typedef std::unordered_map<unsigned int, BranchValVecLayMap> TkToBranchValVecLayMapMap;
typedef TkToBranchValVecLayMapMap::iterator TkToBVVMMIter;
typedef BranchValVecLayMap::iterator BVVLMiter;

// other typedefs
typedef std::unordered_map<unsigned int,const Track*>  TkToTkRefMap;
typedef std::unordered_map<unsigned int,TrackRefVec>   TkToTkRefVecMap;
typedef std::unordered_map<unsigned int,TrackState>    TkToTSMap;   
typedef std::unordered_map<unsigned int,TSVec>         TkToTSVecMap;
typedef std::unordered_map<unsigned int,TSLayerPairVec> TkToTSLayerPairVecMap;

class TTreeValidation : public Validation {
public:
  TTreeValidation(std::string fileName);

  void resetValidationMaps() override;
  
  void collectSimTkTSVecMapInfo(unsigned int mcTrackID, const TSVec& initTSs) override;

  void collectSeedTkCFMapInfo(unsigned int seedID, const TrackState& cfitStateHit0) override;
  void collectSeedTkTSLayerPairVecMapInfo(unsigned int seedID, const TSLayerPairVec& updatedStates) override;

  void collectBranchingInfo(unsigned int seedID, unsigned int ilayer,
                            float nSigmaDeta, float etaBinMinus, unsigned int etaBinPlus,
                            float nSigmaDphi, unsigned int phiBinMinus, unsigned int phiBinPlus,
                            const std::vector<unsigned int>& cand_hit_indices, const std::vector<unsigned int>& cand_hits_branches) override;

  void collectFitTkCFMapInfo(unsigned int seedID, const TrackState& cfitStateHit0) override;
  void collectFitTkTSLayerPairVecMapInfo(unsigned int seedID, const TSLayerPairVec& updatedStates) override;

  void fillSegmentTree(const BinInfoMap& segmentMap, unsigned int evtID) override;

  void fillBranchTree(unsigned int evtID) override;

  void makeSimTkToRecoTksMaps(const TrackVec& evt_seed_tracks, const TrackVec& evt_build_tracks, const TrackVec& evt_fit_tracks) override;
  void mapSimTkToRecoTks(const TrackVec& evt_tracks, TkToTkRefVecMap& simTkMap);
  void fillEffTree(const TrackVec& evt_sim_tracks, unsigned int ievt) override;

  void makeSeedTkToRecoTkMaps(const TrackVec& evt_build_tracks, const TrackVec& evt_fit_tracks) override;
  void mapSeedTkToRecoTk(const TrackVec& evt_tracks, TkToTkRefMap& seedTkMap);
  void fillFakeRateTree(const TrackVec& evt_seed_tracks, unsigned int ievt) override;

  void fillConfigTree(const std::vector<double>& ticks) override;

  void saveTTrees() override;

  TFile* f_;

  // segment map tree
  TTree* segtree_;
  unsigned int evtID_seg_=0,layer_seg_=0,etabin_seg_=0,phibin_seg_=0,nHits_seg_=0;

  // declare vector of trackStates for various track collections to be used in pulls/resolutions
  TkToTSVecMap simTkTSVecMap_;
  
  TkToTSMap seedTkCFMap_;
  TkToTSLayerPairVecMap seedTkTSLayerPairVecMap_;

  TkToBranchValVecLayMapMap seedToBranchValVecLayMapMap_;

  TkToTSMap fitTkCFMap_;
  TkToTSLayerPairVecMap fitTkTSLayerPairVecMap_;

  // build branching tree
  TTree* tree_br_;
  unsigned int evtID_br_=0,seedID_br_=0;
  unsigned int layer_=0,cands_=0;
  unsigned int uniqueEtaPhiBins_,uniqueHits_,uniqueBranches_;
  std::vector<unsigned int> candEtaPhiBins_,candHits_,candBranches_;
  std::vector<float> candnSigmaDeta_,candnSigmaDphi_;
  
  // efficiency trees and variables
  TkToTkRefVecMap simToSeedMap_;
  TkToTkRefVecMap simToBuildMap_;
  TkToTkRefVecMap simToFitMap_;

  TTree* efftree_;  
  unsigned int evtID_eff_=0,mcID_eff_=0;
  int   mcmask_seed_eff_=0,mcmask_build_eff_=0,mcmask_fit_eff_=0;

  int   seedID_seed_eff_=0,seedID_build_eff_=0,seedID_fit_eff_=0;

  // for efficiency and duplicate rate plots
  float pt_mc_gen_eff_=0.,phi_mc_gen_eff_=0.,eta_mc_gen_eff_=0.;
  int   nHits_mc_eff_=0;

  // for position pulls and geo plots
  std::vector<float> x_mc_reco_hit_eff_,y_mc_reco_hit_eff_,z_mc_reco_hit_eff_;
  float x_mc_gen_vrx_eff_=0.,y_mc_gen_vrx_eff_=0.,z_mc_gen_vrx_eff_=0.;

  // for track resolutions / pulls
  float pt_mc_seed_eff_=0.,pt_mc_build_eff_=0.,pt_mc_fit_eff_=0.;
  float pt_seed_eff_=0.,pt_build_eff_=0.,pt_fit_eff_=0.,ept_seed_eff_=0.,ept_build_eff_=0.,ept_fit_eff_=0.;
  float phi_mc_seed_eff_=0.,phi_mc_build_eff_=0.,phi_mc_fit_eff_=0.;
  float phi_seed_eff_=0.,phi_build_eff_=0.,phi_fit_eff_=0.,ephi_seed_eff_=0.,ephi_build_eff_=0.,ephi_fit_eff_=0.;
  float eta_mc_seed_eff_=0.,eta_mc_build_eff_=0.,eta_mc_fit_eff_=0.;
  float eta_seed_eff_=0.,eta_build_eff_=0.,eta_fit_eff_=0.,eeta_seed_eff_=0.,eeta_build_eff_=0.,eeta_fit_eff_=0.;
  
  // for conformal fit resolutions/pulls/residuals
  float x_mc_cf_seed_eff_=0.,y_mc_cf_seed_eff_=0.,z_mc_cf_seed_eff_=0.,x_mc_cf_fit_eff_=0.,y_mc_cf_fit_eff_=0.,z_mc_cf_fit_eff_=0.;
  float x_cf_seed_eff_=0.,y_cf_seed_eff_=0.,z_cf_seed_eff_=0.,x_cf_fit_eff_=0.,y_cf_fit_eff_=0.,z_cf_fit_eff_=0.;
  float ex_cf_seed_eff_=0.,ey_cf_seed_eff_=0.,ez_cf_seed_eff_=0.,ex_cf_fit_eff_=0.,ey_cf_fit_eff_=0.,ez_cf_fit_eff_=0.;

  float px_mc_cf_seed_eff_=0.,py_mc_cf_seed_eff_=0.,pz_mc_cf_seed_eff_=0.,px_mc_cf_fit_eff_=0.,py_mc_cf_fit_eff_=0.,pz_mc_cf_fit_eff_=0.;
  float px_cf_seed_eff_=0.,py_cf_seed_eff_=0.,pz_cf_seed_eff_=0.,px_cf_fit_eff_=0.,py_cf_fit_eff_=0.,pz_cf_fit_eff_=0.;
  float epx_cf_seed_eff_=0.,epy_cf_seed_eff_=0.,epz_cf_seed_eff_=0.,epx_cf_fit_eff_=0.,epy_cf_fit_eff_=0.,epz_cf_fit_eff_=0.;

  float pt_mc_cf_seed_eff_=0.,invpt_mc_cf_seed_eff_=0.,phi_mc_cf_seed_eff_=0.,theta_mc_cf_seed_eff_=0.,pt_mc_cf_fit_eff_=0.,invpt_mc_cf_fit_eff_=0.,phi_mc_cf_fit_eff_=0.,theta_mc_cf_fit_eff_=0.;
  float pt_cf_seed_eff_=0.,invpt_cf_seed_eff_=0.,phi_cf_seed_eff_=0.,theta_cf_seed_eff_=0.,pt_cf_fit_eff_=0.,invpt_cf_fit_eff_=0.,phi_cf_fit_eff_=0.,theta_cf_fit_eff_=0.;
  float ept_cf_seed_eff_=0.,einvpt_cf_seed_eff_=0.,ephi_cf_seed_eff_=0.,etheta_cf_seed_eff_=0.,ept_cf_fit_eff_=0.,einvpt_cf_fit_eff_=0.,ephi_cf_fit_eff_=0.,etheta_cf_fit_eff_=0.;

  // for position pulls
  std::vector<float> x_lay_seed_eff_,y_lay_seed_eff_,z_lay_seed_eff_,x_lay_fit_eff_,y_lay_fit_eff_,z_lay_fit_eff_;
  std::vector<float> ex_lay_seed_eff_,ey_lay_seed_eff_,ez_lay_seed_eff_,ex_lay_fit_eff_,ey_lay_fit_eff_,ez_lay_fit_eff_;
  std::vector<unsigned int> layers_seed_eff_,layers_fit_eff_;

  // for hit countings
  int   nHits_seed_eff_=0,nHits_build_eff_=0,nHits_fit_eff_=0;
  int   nHitsMatched_seed_eff_=0,nHitsMatched_build_eff_=0,nHitsMatched_fit_eff_=0;
  float fracHitsMatched_seed_eff_=0,fracHitsMatched_build_eff_=0,fracHitsMatched_fit_eff_=0;

  float hitchi2_seed_eff_=0.,hitchi2_build_eff_=0.,hitchi2_fit_eff_=0.;
  float helixchi2_seed_eff_=0.,helixchi2_build_eff_=0.,helixchi2_fit_eff_=0.;

  // for duplicate track matches
  int   duplmask_seed_eff_=0,duplmask_build_eff_=0,duplmask_fit_eff_=0;
  int   nTkMatches_seed_eff_=0,nTkMatches_build_eff_=0,nTkMatches_fit_eff_=0;

  // fake rate tree and variables
  TkToTkRefMap seedToBuildMap_;
  TkToTkRefMap seedToFitMap_;
  
  TTree* fakeratetree_;
  unsigned int evtID_FR_=0,seedID_FR_=0;

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

  // configuration tree
  TTree* configtree_;
  float        simtime_=0.,segtime_=0.,seedtime_=0.,buildtime_=0.,fittime_=0.,hlvtime_=0.;
  unsigned int Ntracks_=0,Nevents_=0;
  unsigned int nLayers_=0;
  float        fRadialSpacing_=0.,fRadialExtent_=0.,fInnerSensorSize_=0.,fOuterSensorSize_=0.;
  float        fEtaDet_=0.,fPhiFactor_=0.;
  unsigned int nPhiPart_=0,nEtaPart_=0;
  unsigned int nlayers_per_seed_=0,maxCand_=0;
  float        chi2Cut_=0.,nSigma_=0.,minDPhi_=0.,maxDPhi_=0.,minDEta_=0.,maxDEta_=0.;
  float        beamspotX_=0.,beamspotY_=0.,beamspotZ_=0.;
  float        minSimPt_=0.,maxSimPt_=0.;
  float        hitposerrXY_=0.,hitposerrZ_=0.,hitposerrR_=0.;
  float        varXY_=0.,varZ_=0.;
  unsigned int nTotHit_=0;
  float        ptinverr049_=0.,phierr049_=0.,thetaerr049_=0.,ptinverr012_=0.,phierr012_=0.,thetaerr012_=0.;

  std::mutex glock_;
};
#endif
#endif
