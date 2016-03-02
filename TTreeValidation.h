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

  int etaBinMinus;
  int etaBinPlus;
  int phiBinMinus;
  int phiBinPlus;

  std::vector<int> cand_hit_indices;
  std::vector<int> branch_hit_indices; // size of was just branches
};

// Branch Tree type definitions
typedef std::vector<BranchVal> BranchValVec;
typedef std::unordered_map<int, BranchValVec> BranchValVecLayMap;
typedef std::unordered_map<int, BranchValVecLayMap> TkIDToBranchValVecLayMapMap;
typedef TkIDToBranchValVecLayMapMap::iterator TkIDToBVVMMIter;
typedef BranchValVecLayMap::iterator BVVLMiter;

// Map typedefs needed for mapping different sets of tracks to another
typedef std::unordered_map<int,int>               TkIDToTkIDMap;
typedef std::unordered_map<int,std::vector<int> > TkIDToTkIDVecMap;
typedef std::unordered_map<int,TrackState>        TkIDToTSMap;   
typedef std::unordered_map<int,TSVec>             TkIDToTSVecMap;
typedef std::unordered_map<int,TSLayerPairVec>    TkIDToTSLayerPairVecMap;

class TTreeValidation : public Validation {
public:
  TTreeValidation(std::string fileName);

  void initializeSeedInfoTree();
  void initializeSeedTree();
  void initializeSegmentTree();
  void initializeBranchTree();
  void initializeEfficiencyTree();
  void initializeFakeRateTree();  
  void initializeGeometryTree();
  void initializeConformalTree();
  void initializeConfigTree();
  
  void alignTrackExtra(TrackVec& evt_tracks, TrackExtraVec& evt_extra) override;

  void collectSimTkTSVecMapInfo(int mcTrackID, const TSVec& initTSs) override;
  void collectSeedTkCFMapInfo(int seedID, const TrackState& cfitStateHit0) override;
  void collectSeedTkTSLayerPairVecMapInfo(int seedID, const TSLayerPairVec& updatedStates) override;
  void collectBranchingInfo(int seedID, int ilayer,
                            float nSigmaDeta, float etaBinMinus, int etaBinPlus,
                            float nSigmaDphi, int phiBinMinus, int phiBinPlus,
                            const std::vector<int>& cand_hit_indices, const std::vector<int>& cand_hits_branches) override;
  void collectFitTkCFMapInfo(int seedID, const TrackState& cfitStateHit0) override;
  void collectFitTkTSLayerPairVecMapInfo(int seedID, const TSLayerPairVec& updatedStates) override;

  void resetValidationMaps() override;
  void makeSimTkToRecoTksMaps(Event& ev) override;
  void mapSimTkToRecoTks(const TrackVec& evt_tracks, TrackExtraVec& evt_extra, const std::vector<HitVec>& layerHits, 
			 const MCHitInfoVec&, TkIDToTkIDVecMap& simTkMap);
  void makeSeedTkToRecoTkMaps(Event& ev) override;
  void mapSeedTkToRecoTk(const TrackVec& evt_tracks, const TrackExtraVec& evt_extras, TkIDToTkIDMap& seedTkMap);

  void fillSeedInfoTree(const TripletIdxVec& hit_triplets, const Event& ev) override;
  void fillSeedTree(const TripletIdxVec& hit_triplets, const TripletIdxVec& filtered_triplets, const Event& ev) override;
  
  void fillSegmentTree(const BinInfoMap& segmentMap, int evtID) override;
  void fillBranchTree(int evtID) override;
  void fillEfficiencyTree(const Event& ev) override;
  void fillFakeRateTree(const Event& ev) override;
  void fillGeometryTree(const Event& ev) override;
  void fillConformalTree(const Event& ev) override;
  void fillConfigTree(const std::vector<double>& ticks) override;

  void saveTTrees() override;

 private:
  TFile* f_; // output file!

  TkIDToTSLayerPairVecMap seedTkTSLayerPairVecMap_; // used for position pulls for seed track
  TkIDToTSLayerPairVecMap fitTkTSLayerPairVecMap_; //  used for position pulls for fit track
  TkIDToBranchValVecLayMapMap seedToBranchValVecLayMapMap_; // map created inside collectBranchingInfo

  TkIDToTSVecMap simTkTSVecMap_; // used for pulls (map all sim track TS to sim ID)
  TkIDToTSMap seedTkCFMap_; // map CF TS to seedID of seed track
  TkIDToTSMap fitTkCFMap_; // map CF TS to seedID of fit track

  // Sim to Reco Maps
  TkIDToTkIDVecMap simToSeedMap_;
  TkIDToTkIDVecMap simToBuildMap_;
  TkIDToTkIDVecMap simToFitMap_;

  // Reco to Reco Maps
  TkIDToTkIDMap seedToBuildMap_;
  TkIDToTkIDMap seedToFitMap_;

  // seedinfo tree
  TTree* seedinfotree_;
  int evtID_seedinfo_;
  int mcID_seedinfo_;
  float pt_gen_seedinfo_;
  float eta_gen_seedinfo_;
  float phi_gen_seedinfo_;
  float a,b,r,d0;
  bool pass;

  // seed tree
  TTree* seedtree_;
  int evtID_seed_;
  int nTkAll_;
  int nTkAllMC_;
  int nTkCut_;
  int nTkCutMC_;

  // segment map tree
  TTree* segtree_;
  int   evtID_seg_=0,layer_seg_=0,etabin_seg_=0,phibin_seg_=0,nHits_seg_=0;

  // build branching tree
  TTree* tree_br_;
  int   evtID_br_=0,seedID_br_=0;
  int   layer_=0,cands_=0;
  int   uniqueEtaPhiBins_,uniqueHits_,uniqueBranches_;
  std::vector<int> candEtaPhiBins_,candHits_,candBranches_;
  std::vector<float> candnSigmaDeta_,candnSigmaDphi_;
  
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

  // Geometry Tree 
  TTree* geotree_;  
  int   evtID_geo_=0,mcID_geo_=0;
  int   mcmask_seed_geo_=0,mcmask_build_geo_=0,mcmask_fit_geo_=0;
  int   seedID_seed_geo_=0,seedID_build_geo_=0,seedID_fit_geo_=0;

  float x_mc_gen_vrx_geo_=0.,y_mc_gen_vrx_geo_=0.,z_mc_gen_vrx_geo_=0.;
  std::vector<int>   layers_seed_geo_,layers_fit_geo_;
  std::vector<float> x_mc_reco_hit_geo_,y_mc_reco_hit_geo_,z_mc_reco_hit_geo_;
  std::vector<float> x_lay_seed_geo_,y_lay_seed_geo_,z_lay_seed_geo_,x_lay_fit_geo_,y_lay_fit_geo_,z_lay_fit_geo_;
  std::vector<float> ex_lay_seed_geo_,ey_lay_seed_geo_,ez_lay_seed_geo_,ex_lay_fit_geo_,ey_lay_fit_geo_,ez_lay_fit_geo_;

  // Conformal Tree 
  TTree* cftree_;  
  int   evtID_cf_=0,mcID_cf_=0;
  int   mcmask_seed_cf_=0,mcmask_build_cf_=0,mcmask_fit_cf_=0;
  int   seedID_seed_cf_=0,seedID_build_cf_=0,seedID_fit_cf_=0;

  float x_mc_seed_cf_=0.,y_mc_seed_cf_=0.,z_mc_seed_cf_=0.,x_mc_fit_cf_=0.,y_mc_fit_cf_=0.,z_mc_fit_cf_=0.;
  float x_seed_cf_=0.,y_seed_cf_=0.,z_seed_cf_=0.,x_fit_cf_=0.,y_fit_cf_=0.,z_fit_cf_=0.;
  float ex_seed_cf_=0.,ey_seed_cf_=0.,ez_seed_cf_=0.,ex_fit_cf_=0.,ey_fit_cf_=0.,ez_fit_cf_=0.;

  float px_mc_seed_cf_=0.,py_mc_seed_cf_=0.,pz_mc_seed_cf_=0.,px_mc_fit_cf_=0.,py_mc_fit_cf_=0.,pz_mc_fit_cf_=0.;
  float px_seed_cf_=0.,py_seed_cf_=0.,pz_seed_cf_=0.,px_fit_cf_=0.,py_fit_cf_=0.,pz_fit_cf_=0.;
  float epx_seed_cf_=0.,epy_seed_cf_=0.,epz_seed_cf_=0.,epx_fit_cf_=0.,epy_fit_cf_=0.,epz_fit_cf_=0.;

  float pt_mc_seed_cf_=0.,invpt_mc_seed_cf_=0.,phi_mc_seed_cf_=0.,theta_mc_seed_cf_=0.,pt_mc_fit_cf_=0.,
        invpt_mc_fit_cf_=0.,phi_mc_fit_cf_=0.,theta_mc_fit_cf_=0.;
  float pt_seed_cf_=0.,invpt_seed_cf_=0.,phi_seed_cf_=0.,theta_seed_cf_=0.,pt_fit_cf_=0.,invpt_fit_cf_=0.,phi_fit_cf_=0.,theta_fit_cf_=0.;
  float ept_seed_cf_=0.,einvpt_seed_cf_=0.,ephi_seed_cf_=0.,etheta_seed_cf_=0.,ept_fit_cf_=0.,einvpt_fit_cf_=0.,ephi_fit_cf_=0.,etheta_fit_cf_=0.;

  // Configuration tree
  TTree* configtree_;
  float simtime_=0.,segtime_=0.,seedtime_=0.,buildtime_=0.,fittime_=0.,hlvtime_=0.;
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

  std::mutex glock_;
};
#endif
#endif
