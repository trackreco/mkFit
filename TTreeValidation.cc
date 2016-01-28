// N.B. Mask assignments
// --> mcmask_[reco] == 1, "associated" reco to sim track [possible duplmask_[reco] == 1,0] {eff and FR}
// --> mcmask_[reco] == 0, "unassociated" reco to sim track. by definition no duplicates (no reco to associate to sim tracks!) [possible duplmask_[reco] == 2 {eff and FR}]
// --> mcmask_[reco] == -1, "no matching seed to build/fit" track, therefore no build/fit track to match sim! [possible duplmask_[reco] == -1] {FR only} 

// --> nTkMatches_[reco] > 1,   n reco tracks associated to the same sim track ID {eff only}
// --> nTkMatches_[reco] == 1,  1 reco track associated to single sim track ID {eff only}
// --> nTkMatches_[reco] == -99, no reco to sim match {eff only}

// excluding position variables, as position could be -99!
// --> reco var == -99, "unassociated" reco to sim track, mcTrackID == 999999 [possible mcmask_[reco] == 0; possible duplmask_[reco] == 2] {eff only}
// --> sim  var == -99, "unassociated" reco to sim track, mcTrackID == 999999 [possible mcmask_[reco] == 0; possible duplmask_[reco] == 2] {FR only}
// --> reco/sim var == -100, "no matching seed to build/fit" track, fill all reco/sim variables -100 [possible mcmask_[reco] == -1, possible duplmask_[reco] == -1] {FR only}

// --> seedmask_[reco] == 1, matching seed to reco/fit track [possible mcmask_[reco] == 0,1; possible duplmask_[reco] == 0,1,2] {FR only}
// --> seedmask_[reco] == 0, no matching seed to reco/fit track [possible mcmask_[reco] == -1; possible duplmask_[reco] == -1] {FR only}

// --> duplmask_[reco] == 0, only "associated" reco to sim track [possible mcmask_[reco] == 1] {eff and FR}
// --> duplmask_[reco] == 1, more than one "associated" reco to sim track [possible mcmask_[reco] == 1] {eff and FR}
// --> duplmask_[reco] == 2, no "associated" reco to sim track [possible mcmask_[reco] == 0] {eff and FR}
// --> duplmask_[reco] == -1, no matching built/fit track for given seed [possible mcmask_[reco] == -1] {FR only}

// --> reco var == -10, variable not yet implemented for given track object

// position reco variables
// --> layers_[reco]    ==  -1, reco unassociated to sim tk {eff only}
// --> reco pos+err var == -2000, reco tk is unassociated to sim tk, just filled once {eff only}

#include "TTreeValidation.h"
#include "Event.h"
#include "Config.h"
#include "Propagation.h"
#ifndef NO_ROOT

inline bool sortByHitsChi2(const Track & cand1, const Track & cand2)
{
  if (cand1.nFoundHits()==cand2.nFoundHits()) return cand1.chi2()<cand2.chi2();
  return cand1.nFoundHits()>cand2.nFoundHits();
}

inline float computeHelixChi2(const SVector6& simParams, const SVector6& recoParams, const SMatrixSym66& recoErrs)
{ 
  float chi2 = 0;
  for(auto i = 0; i < Config::nParams; i++){
    float delta = simParams.At(i) - recoParams.At(i);
    chi2 += (delta*delta) / recoErrs.At(i,i);
  }
  return chi2 / (Config::nParams - 1);
}

TTreeValidation::TTreeValidation(std::string fileName)
{
  std::lock_guard<std::mutex> locker(glock_);
  gROOT->ProcessLine("#include <vector>");
  f_ = TFile::Open(fileName.c_str(), "recreate");

  initializeSegmentTree();
  initializeBranchTree();
  initializeEfficiencyTree();
  initializeFakeRateTree();  
  initializeGeometryTree();
  initializeConformalTree();
  initializeConfigTree();
}

void TTreeValidation::initializeSegmentTree(){
  // segment validation
  segtree_ = new TTree("segtree","segtree");
  segtree_->Branch("evtID",&evtID_seg_);
  segtree_->Branch("layer",&layer_seg_);
  segtree_->Branch("etabin",&etabin_seg_);
  segtree_->Branch("phibin",&phibin_seg_);
  segtree_->Branch("nHits",&nHits_seg_);
}

void TTreeValidation::initializeBranchTree(){
  // build validation
  tree_br_ = new TTree("tree_br","tree_br");
  tree_br_->Branch("evtID",&evtID_br_);
  tree_br_->Branch("seedID",&seedID_br_);

  tree_br_->Branch("layer",&layer_);
  tree_br_->Branch("cands",&cands_);

  tree_br_->Branch("candEtaPhiBins",&candEtaPhiBins_);
  tree_br_->Branch("candHits",&candHits_);
  tree_br_->Branch("candBranches",&candBranches_);
  
  tree_br_->Branch("uniqueEtaPhiBins",&uniqueEtaPhiBins_);
  tree_br_->Branch("uniqueHits",&uniqueHits_);
  tree_br_->Branch("uniqueBranches",&uniqueBranches_);
      
  tree_br_->Branch("candnSigmaDeta",&candnSigmaDeta_);
  tree_br_->Branch("candnSigmaDphi",&candnSigmaDphi_);
}

void TTreeValidation::initializeEfficiencyTree(){  
  // efficiency validation
  efftree_ = new TTree("efftree","efftree");
  efftree_->Branch("evtID",&evtID_eff_);
  efftree_->Branch("mcID",&mcID_eff_);

  efftree_->Branch("nHits_mc",&nHits_mc_eff_);

  efftree_->Branch("seedID_seed",&seedID_seed_eff_);
  efftree_->Branch("seedID_build",&seedID_build_eff_);
  efftree_->Branch("seedID_fit",&seedID_fit_eff_);

  efftree_->Branch("pt_mc_gen",&pt_mc_gen_eff_);
  efftree_->Branch("phi_mc_gen",&phi_mc_gen_eff_);
  efftree_->Branch("eta_mc_gen",&eta_mc_gen_eff_);

  efftree_->Branch("mcmask_seed",&mcmask_seed_eff_);
  efftree_->Branch("mcmask_build",&mcmask_build_eff_);
  efftree_->Branch("mcmask_fit",&mcmask_fit_eff_);

  efftree_->Branch("pt_mc_seed",&pt_mc_seed_eff_);
  efftree_->Branch("pt_seed",&pt_seed_eff_);
  efftree_->Branch("ept_seed",&ept_seed_eff_);
  efftree_->Branch("pt_mc_build",&pt_mc_build_eff_);
  efftree_->Branch("pt_build",&pt_build_eff_);
  efftree_->Branch("ept_build",&ept_build_eff_);
  efftree_->Branch("pt_mc_fit",&pt_mc_fit_eff_);
  efftree_->Branch("pt_fit",&pt_fit_eff_);
  efftree_->Branch("ept_fit",&ept_fit_eff_);

  efftree_->Branch("phi_mc_seed",&phi_mc_seed_eff_);
  efftree_->Branch("phi_seed",&phi_seed_eff_);
  efftree_->Branch("ephi_seed",&ephi_seed_eff_);
  efftree_->Branch("phi_mc_build",&phi_mc_build_eff_);
  efftree_->Branch("phi_build",&phi_build_eff_);
  efftree_->Branch("ephi_build",&ephi_build_eff_);
  efftree_->Branch("phi_mc_fit",&phi_mc_fit_eff_);
  efftree_->Branch("phi_fit",&phi_fit_eff_);
  efftree_->Branch("ephi_fit",&ephi_fit_eff_);

  efftree_->Branch("eta_mc_seed",&eta_mc_seed_eff_);
  efftree_->Branch("eta_seed",&eta_seed_eff_);
  efftree_->Branch("eeta_seed",&eeta_seed_eff_);
  efftree_->Branch("eta_mc_build",&eta_mc_build_eff_);
  efftree_->Branch("eta_build",&eta_build_eff_);
  efftree_->Branch("eeta_build",&eeta_build_eff_);
  efftree_->Branch("eta_mc_fit",&eta_mc_fit_eff_);
  efftree_->Branch("eta_fit",&eta_fit_eff_);
  efftree_->Branch("eeta_fit",&eeta_fit_eff_);

  efftree_->Branch("nHits_seed",&nHits_seed_eff_);
  efftree_->Branch("nHits_build",&nHits_build_eff_);
  efftree_->Branch("nHits_fit",&nHits_fit_eff_);

  efftree_->Branch("nHitsMatched_seed",&nHitsMatched_seed_eff_);
  efftree_->Branch("nHitsMatched_build",&nHitsMatched_build_eff_);
  efftree_->Branch("nHitsMatched_fit",&nHitsMatched_fit_eff_);

  efftree_->Branch("fracHitsMatched_seed",&fracHitsMatched_seed_eff_);
  efftree_->Branch("fracHitsMatched_build",&fracHitsMatched_build_eff_);
  efftree_->Branch("fracHitsMatched_fit",&fracHitsMatched_fit_eff_);

  efftree_->Branch("hitchi2_seed",&hitchi2_seed_eff_);
  efftree_->Branch("hitchi2_build",&hitchi2_build_eff_);
  efftree_->Branch("hitchi2_fit",&hitchi2_fit_eff_);

  efftree_->Branch("helixchi2_seed",&helixchi2_seed_eff_);
  efftree_->Branch("helixchi2_build",&helixchi2_build_eff_);
  efftree_->Branch("helixchi2_fit",&helixchi2_fit_eff_);

  efftree_->Branch("duplmask_seed",&duplmask_seed_eff_);
  efftree_->Branch("duplmask_build",&duplmask_build_eff_);
  efftree_->Branch("duplmask_fit",&duplmask_fit_eff_);

  efftree_->Branch("nTkMatches_seed",&nTkMatches_seed_eff_);
  efftree_->Branch("nTkMatches_build",&nTkMatches_build_eff_);
  efftree_->Branch("nTkMatches_fit",&nTkMatches_fit_eff_);
}

void TTreeValidation::initializeFakeRateTree(){
  // fake rate validation
  fakeratetree_ = new TTree("fakeratetree","fakeratetree");

  fakeratetree_->Branch("evtID",&evtID_FR_);
  fakeratetree_->Branch("seedID",&seedID_FR_);

  fakeratetree_->Branch("seedmask_seed",&seedmask_seed_FR_);
  fakeratetree_->Branch("seedmask_build",&seedmask_build_FR_);
  fakeratetree_->Branch("seedmask_fit",&seedmask_fit_FR_);

  fakeratetree_->Branch("pt_seed",&pt_seed_FR_);
  fakeratetree_->Branch("ept_seed",&ept_seed_FR_);
  fakeratetree_->Branch("pt_build",&pt_build_FR_);
  fakeratetree_->Branch("ept_build",&ept_build_FR_);
  fakeratetree_->Branch("pt_fit",&pt_fit_FR_);
  fakeratetree_->Branch("ept_fit",&ept_fit_FR_);

  fakeratetree_->Branch("phi_seed",&phi_seed_FR_);
  fakeratetree_->Branch("ephi_seed",&ephi_seed_FR_);
  fakeratetree_->Branch("phi_build",&phi_build_FR_);
  fakeratetree_->Branch("ephi_build",&ephi_build_FR_);
  fakeratetree_->Branch("phi_fit",&phi_fit_FR_);
  fakeratetree_->Branch("ephi_fit",&ephi_fit_FR_);

  fakeratetree_->Branch("eta_seed",&eta_seed_FR_);
  fakeratetree_->Branch("eeta_seed",&eeta_seed_FR_);
  fakeratetree_->Branch("eta_build",&eta_build_FR_);
  fakeratetree_->Branch("eeta_build",&eeta_build_FR_);
  fakeratetree_->Branch("eta_fit",&eta_fit_FR_);
  fakeratetree_->Branch("eeta_fit",&eeta_fit_FR_);

  fakeratetree_->Branch("nHits_seed",&nHits_seed_FR_);
  fakeratetree_->Branch("nHits_build",&nHits_build_FR_);
  fakeratetree_->Branch("nHits_fit",&nHits_fit_FR_);

  fakeratetree_->Branch("nHitsMatched_seed",&nHitsMatched_seed_FR_);
  fakeratetree_->Branch("nHitsMatched_build",&nHitsMatched_build_FR_);
  fakeratetree_->Branch("nHitsMatched_fit",&nHitsMatched_fit_FR_);

  fakeratetree_->Branch("fracHitsMatched_seed",&fracHitsMatched_seed_FR_);
  fakeratetree_->Branch("fracHitsMatched_build",&fracHitsMatched_build_FR_);
  fakeratetree_->Branch("fracHitsMatched_fit",&fracHitsMatched_fit_FR_);

  fakeratetree_->Branch("hitchi2_seed",&hitchi2_seed_FR_);
  fakeratetree_->Branch("hitchi2_build",&hitchi2_build_FR_);
  fakeratetree_->Branch("hitchi2_fit",&hitchi2_fit_FR_);

  // sim info of seed,build,fit tracks
  fakeratetree_->Branch("mcID_seed",&mcID_seed_FR_);
  fakeratetree_->Branch("mcID_build",&mcID_build_FR_);
  fakeratetree_->Branch("mcID_fit",&mcID_fit_FR_);
  
  fakeratetree_->Branch("mcmask_seed",&mcmask_seed_FR_);
  fakeratetree_->Branch("mcmask_build",&mcmask_build_FR_);
  fakeratetree_->Branch("mcmask_fit",&mcmask_fit_FR_);

  fakeratetree_->Branch("pt_mc_seed",&pt_mc_seed_FR_);
  fakeratetree_->Branch("pt_mc_build",&pt_mc_build_FR_);
  fakeratetree_->Branch("pt_mc_fit",&pt_mc_fit_FR_);

  fakeratetree_->Branch("phi_mc_seed",&phi_mc_seed_FR_);
  fakeratetree_->Branch("phi_mc_build",&phi_mc_build_FR_);
  fakeratetree_->Branch("phi_mc_fit",&phi_mc_fit_FR_);

  fakeratetree_->Branch("eta_mc_seed",&eta_mc_seed_FR_);
  fakeratetree_->Branch("eta_mc_build",&eta_mc_build_FR_);
  fakeratetree_->Branch("eta_mc_fit",&eta_mc_fit_FR_);

  fakeratetree_->Branch("nHits_mc_seed",&nHits_mc_seed_FR_);
  fakeratetree_->Branch("nHits_mc_build",&nHits_mc_build_FR_);
  fakeratetree_->Branch("nHits_mc_fit",&nHits_mc_fit_FR_);

  fakeratetree_->Branch("helixchi2_seed",&helixchi2_seed_FR_);
  fakeratetree_->Branch("helixchi2_build",&helixchi2_build_FR_);
  fakeratetree_->Branch("helixchi2_fit",&helixchi2_fit_FR_);

  fakeratetree_->Branch("duplmask_seed",&duplmask_seed_FR_);
  fakeratetree_->Branch("duplmask_build",&duplmask_build_FR_);
  fakeratetree_->Branch("duplmask_fit",&duplmask_fit_FR_);

  fakeratetree_->Branch("iTkMatches_seed",&iTkMatches_seed_FR_);
  fakeratetree_->Branch("iTkMatches_build",&iTkMatches_build_FR_);
  fakeratetree_->Branch("iTkMatches_fit",&iTkMatches_fit_FR_);
}

void TTreeValidation::initializeGeometryTree(){
  // Geometry validation
  geotree_ = new TTree("geotree","geotree");

  geotree_->Branch("evtID",&evtID_geo_);
  geotree_->Branch("mcID",&mcID_geo_);
  geotree_->Branch("seedID_seed",&seedID_seed_geo_);
  geotree_->Branch("seedID_fit",&seedID_fit_geo_);
  geotree_->Branch("mcmask_seed",&mcmask_seed_geo_);
  geotree_->Branch("mcmask_fit",&mcmask_fit_geo_);

  // Generated vertex of MC tracks  
  geotree_->Branch("x_mc_gen_vrx",&x_mc_gen_vrx_geo_);
  geotree_->Branch("y_mc_gen_vrx",&y_mc_gen_vrx_geo_);
  geotree_->Branch("z_mc_gen_vrx",&z_mc_gen_vrx_geo_);
  
  // Full reco hits (vector) from simulation
  geotree_->Branch("x_mc_reco_hit",&x_mc_reco_hit_geo_);
  geotree_->Branch("y_mc_reco_hit",&y_mc_reco_hit_geo_);
  geotree_->Branch("z_mc_reco_hit",&z_mc_reco_hit_geo_);
  
  // Position pull info from reco tracks
  geotree_->Branch("layers_seed",&layers_seed_geo_);
  geotree_->Branch("layers_fit",&layers_fit_geo_);

  geotree_->Branch("x_lay_seed",&x_lay_seed_geo_);
  geotree_->Branch("y_lay_seed",&y_lay_seed_geo_);
  geotree_->Branch("z_lay_seed",&z_lay_seed_geo_);
  geotree_->Branch("ex_lay_seed",&ex_lay_seed_geo_);
  geotree_->Branch("ey_lay_seed",&ey_lay_seed_geo_);
  geotree_->Branch("ez_lay_seed",&ez_lay_seed_geo_);

  geotree_->Branch("x_lay_fit",&x_lay_fit_geo_);
  geotree_->Branch("y_lay_fit",&y_lay_fit_geo_);
  geotree_->Branch("z_lay_fit",&z_lay_fit_geo_);
  geotree_->Branch("ex_lay_fit",&ex_lay_fit_geo_);
  geotree_->Branch("ey_lay_fit",&ey_lay_fit_geo_);
  geotree_->Branch("ez_lay_fit",&ez_lay_fit_geo_);
}

void TTreeValidation::initializeConformalTree(){
  // Conformal Fit validation
  cftree_ = new TTree("cftree","cftree");

  cftree_->Branch("evtID",&evtID_cf_);
  cftree_->Branch("mcID",&mcID_cf_);
  cftree_->Branch("seedID_seed",&seedID_seed_cf_);
  cftree_->Branch("seedID_fit",&seedID_fit_cf_);
  cftree_->Branch("mcmask_seed",&mcmask_seed_cf_);
  cftree_->Branch("mcmask_fit",&mcmask_fit_cf_);

  cftree_->Branch("x_mc_cf_seed",&x_mc_seed_cf_);
  cftree_->Branch("x_cf_seed",&x_seed_cf_);
  cftree_->Branch("ex_cf_seed",&ex_seed_cf_);
  cftree_->Branch("x_mc_cf_fit",&x_mc_fit_cf_);
  cftree_->Branch("x_cf_fit",&x_fit_cf_);
  cftree_->Branch("ex_cf_fit",&ex_fit_cf_);

  cftree_->Branch("y_mc_cf_seed",&y_mc_seed_cf_);
  cftree_->Branch("y_cf_seed",&y_seed_cf_);
  cftree_->Branch("ey_cf_seed",&ey_seed_cf_);
  cftree_->Branch("y_mc_cf_fit",&y_mc_fit_cf_);
  cftree_->Branch("y_cf_fit",&y_fit_cf_);
  cftree_->Branch("ey_cf_fit",&ey_fit_cf_);

  cftree_->Branch("z_mc_cf_seed",&z_mc_seed_cf_);
  cftree_->Branch("z_cf_seed",&z_seed_cf_);
  cftree_->Branch("ez_cf_seed",&ez_seed_cf_);
  cftree_->Branch("z_mc_cf_fit",&z_mc_fit_cf_);
  cftree_->Branch("z_cf_fit",&z_fit_cf_);
  cftree_->Branch("ez_cf_fit",&ez_fit_cf_);

  cftree_->Branch("px_mc_cf_seed",&px_mc_seed_cf_);
  cftree_->Branch("px_cf_seed",&px_seed_cf_);
  cftree_->Branch("epx_cf_seed",&epx_seed_cf_);
  cftree_->Branch("px_mc_cf_fit",&px_mc_fit_cf_);
  cftree_->Branch("px_cf_fit",&px_fit_cf_);
  cftree_->Branch("epx_cf_fit",&epx_fit_cf_);

  cftree_->Branch("py_mc_cf_seed",&py_mc_seed_cf_);
  cftree_->Branch("py_cf_seed",&py_seed_cf_);
  cftree_->Branch("epy_cf_seed",&epy_seed_cf_);
  cftree_->Branch("py_mc_cf_fit",&py_mc_fit_cf_);
  cftree_->Branch("py_cf_fit",&py_fit_cf_);
  cftree_->Branch("epy_cf_fit",&epy_fit_cf_);

  cftree_->Branch("pz_mc_cf_seed",&pz_mc_seed_cf_);
  cftree_->Branch("pz_cf_seed",&pz_seed_cf_);
  cftree_->Branch("epz_cf_seed",&epz_seed_cf_);
  cftree_->Branch("pz_mc_cf_fit",&pz_mc_fit_cf_);
  cftree_->Branch("pz_cf_fit",&pz_fit_cf_);
  cftree_->Branch("epz_cf_fit",&epz_fit_cf_);

  cftree_->Branch("pt_mc_cf_seed",&pt_mc_seed_cf_);
  cftree_->Branch("pt_cf_seed",&pt_seed_cf_);
  cftree_->Branch("ept_cf_seed",&ept_seed_cf_);
  cftree_->Branch("pt_mc_cf_fit",&pt_mc_fit_cf_);
  cftree_->Branch("pt_cf_fit",&pt_fit_cf_);
  cftree_->Branch("ept_cf_fit",&ept_fit_cf_);

  cftree_->Branch("invpt_mc_cf_seed",&invpt_mc_seed_cf_);
  cftree_->Branch("invpt_cf_seed",&invpt_seed_cf_);
  cftree_->Branch("einvpt_cf_seed",&einvpt_seed_cf_);
  cftree_->Branch("invpt_mc_cf_fit",&invpt_mc_fit_cf_);
  cftree_->Branch("invpt_cf_fit",&invpt_fit_cf_);
  cftree_->Branch("einvpt_cf_fit",&einvpt_fit_cf_);

  cftree_->Branch("phi_mc_cf_seed",&phi_mc_seed_cf_);
  cftree_->Branch("phi_cf_seed",&phi_seed_cf_);
  cftree_->Branch("ephi_cf_seed",&ephi_seed_cf_);
  cftree_->Branch("phi_mc_cf_fit",&phi_mc_fit_cf_);
  cftree_->Branch("phi_cf_fit",&phi_fit_cf_);
  cftree_->Branch("ephi_cf_fit",&ephi_fit_cf_);
  
  cftree_->Branch("theta_mc_cf_seed",&theta_mc_seed_cf_);
  cftree_->Branch("theta_cf_seed",&theta_seed_cf_);
  cftree_->Branch("etheta_cf_seed",&etheta_seed_cf_);
  cftree_->Branch("theta_mc_cf_fit",&theta_mc_fit_cf_);
  cftree_->Branch("theta_cf_fit",&theta_fit_cf_);
  cftree_->Branch("etheta_cf_fit",&etheta_fit_cf_);
}

void TTreeValidation::initializeConfigTree(){
  // include config ++ real seeding parameters ...
  configtree_ = new TTree("configtree","configtree");
  configtree_->Branch("simtime",&simtime_);
  configtree_->Branch("segtime",&segtime_);
  configtree_->Branch("seedtime",&seedtime_);
  configtree_->Branch("buildtime",&buildtime_);
  configtree_->Branch("fittime",&fittime_);
  configtree_->Branch("hlvtime",&hlvtime_);

  configtree_->Branch("Ntracks",&Ntracks_);
  configtree_->Branch("Nevents",&Nevents_);
  
  configtree_->Branch("nLayers",&nLayers_);
  configtree_->Branch("fRadialSpacing",&fRadialSpacing_);
  configtree_->Branch("fRadialExtent",&fRadialExtent_);
  configtree_->Branch("fInnerSensorSize",&fInnerSensorSize_);
  configtree_->Branch("fOuterSensorSize",&fOuterSensorSize_);
  configtree_->Branch("fEtaDet",&fEtaDet_);

  configtree_->Branch("nPhiPart",&nPhiPart_);
  configtree_->Branch("fPhiFactor",&fPhiFactor_);
  configtree_->Branch("nEtaPart",&nEtaPart_);

  configtree_->Branch("nlayers_per_seed",&nlayers_per_seed_);
  configtree_->Branch("maxCand",&maxCand_);
  configtree_->Branch("chi2Cut",&chi2Cut_);
  configtree_->Branch("nSigma",&nSigma_);
  configtree_->Branch("minDPhi",&minDPhi_);
  configtree_->Branch("maxDPhi",&maxDPhi_);
  configtree_->Branch("minDEta",&minDEta_);
  configtree_->Branch("maxDEta",&maxDEta_);

  configtree_->Branch("beamspotX",&beamspotX_);
  configtree_->Branch("beamspotY",&beamspotY_);
  configtree_->Branch("beamspotZ",&beamspotZ_);

  configtree_->Branch("minSimPt",&minSimPt_);
  configtree_->Branch("maxSimPt",&maxSimPt_);

  configtree_->Branch("hitposerrXY",&hitposerrXY_);
  configtree_->Branch("hitposerrZ",&hitposerrZ_);
  configtree_->Branch("hitposerrR",&hitposerrR_);

  configtree_->Branch("varXY",&varXY_);
  configtree_->Branch("varZ",&varZ_);

  configtree_->Branch("nTotHit",&nTotHit_);

  configtree_->Branch("ptinverr049",&ptinverr049_);
  configtree_->Branch("phierr049",&phierr049_);
  configtree_->Branch("thetaerr049",&thetaerr049_);
  configtree_->Branch("ptinverr012",&ptinverr012_);
  configtree_->Branch("phierr012",&phierr012_);
  configtree_->Branch("thetaerr012",&thetaerr012_);
}

void TTreeValidation::alignTrackExtra(TrackVec& evt_tracks, TrackExtraVec& evt_extras){
  TrackExtraVec trackExtra_tmp;

  // align temporary tkExVec with new track collection ordering
  for (auto&& track : evt_tracks){ 
    trackExtra_tmp.push_back(evt_extras[track.label()]); // label is old seedID!
  }

  // now copy the temporary back in the old one
  evt_extras = trackExtra_tmp;
  
  // redo track labels to match index in vector
  for (int i = 0; i < evt_tracks.size(); i++){
    evt_tracks[i].setLabel(i);
  }
}

void TTreeValidation::collectSimTkTSVecMapInfo(int mcTrackID, const TSVec& initTSs){
  simTkTSVecMap_[mcTrackID] = initTSs;
}

void TTreeValidation::collectSeedTkCFMapInfo(int seedID, const TrackState& cfitStateHit0){
  seedTkCFMap_[seedID] = cfitStateHit0;
}

void TTreeValidation::collectSeedTkTSLayerPairVecMapInfo(int seedID, const TSLayerPairVec& updatedStates){
  seedTkTSLayerPairVecMap_[seedID] = updatedStates;
}

void TTreeValidation::collectBranchingInfo(int seedID, int ilayer, 
					   float nSigmaDeta, float etaBinMinus, int etaBinPlus, 
					   float nSigmaDphi, int phiBinMinus, int phiBinPlus, 
					   const std::vector<int>& cand_hit_indices, const std::vector<int>& branch_hit_indices){
  
  BranchVal tmpBranchVal;
  tmpBranchVal.nSigmaDeta  = nSigmaDeta;
  tmpBranchVal.etaBinMinus = etaBinMinus;
  tmpBranchVal.etaBinPlus  = etaBinPlus;
  tmpBranchVal.nSigmaDphi  = nSigmaDphi;
  tmpBranchVal.phiBinMinus = phiBinMinus;
  tmpBranchVal.phiBinPlus  = phiBinPlus;
  tmpBranchVal.cand_hit_indices   = cand_hit_indices;
  tmpBranchVal.branch_hit_indices = branch_hit_indices; // size of vector formerly known as just "branches"
  
  seedToBranchValVecLayMapMap_[seedID][ilayer].push_back(tmpBranchVal);
}

void TTreeValidation::collectFitTkCFMapInfo(int seedID, const TrackState& cfitStateHit0){
  fitTkCFMap_[seedID] = cfitStateHit0;
}

void TTreeValidation::collectFitTkTSLayerPairVecMapInfo(int seedID, const TSLayerPairVec& updatedStates){
  fitTkTSLayerPairVecMap_[seedID] = updatedStates;
}

void TTreeValidation::resetValidationMaps(){
  std::lock_guard<std::mutex> locker(glock_);
  
  // reset map of sim tracks to MC track states (pre-smearing)
  simTkTSVecMap_.clear(); 

  // reset map of reco tracks to Conformal Fit track states 
  seedTkCFMap_.clear();
  fitTkCFMap_.clear(); 

  // reset map of reco tracks to updated layer track states and bool pairs
  seedTkTSLayerPairVecMap_.clear();
  fitTkTSLayerPairVecMap_.clear();

  // reset branching map
  seedToBranchValVecLayMapMap_.clear();

  // reset map of sim tracks to reco tracks
  simToSeedMap_.clear();
  simToBuildMap_.clear();
  simToFitMap_.clear();

  // reset map of seed tracks to reco tracks
  seedToBuildMap_.clear();
  seedToFitMap_.clear();
}

void TTreeValidation::makeSimTkToRecoTksMaps(Event& ev){
  std::lock_guard<std::mutex> locker(glock_);
  // set mcTkIDs... and sort by each (simTracks set in order by default!)
  mapSimTkToRecoTks(ev.seedTracks_,ev.seedTracksExtra_,ev.layerHits_,ev.simHitsInfo_,simToSeedMap_);
  mapSimTkToRecoTks(ev.candidateTracks_,ev.candidateTracksExtra_,ev.layerHits_,ev.simHitsInfo_,simToBuildMap_);
  mapSimTkToRecoTks(ev.fitTracks_,ev.fitTracksExtra_,ev.layerHits_,ev.simHitsInfo_,simToFitMap_);
}

void TTreeValidation::mapSimTkToRecoTks(const TrackVec& evt_tracks, TrackExtraVec& evt_extras, const std::vector<HitVec>& layerHits, 
					const MCHitInfoVec& mcHitInfo, TkIDToTkIDVecMap& simTkMap){
  for (auto itrack = 0; itrack < evt_tracks.size(); ++itrack){
    auto&& track(evt_tracks[itrack]);
    auto&& extra(evt_extras[itrack]);
    extra.setMCTrackIDInfo(track, layerHits, mcHitInfo);
    if (extra.mcTrackID() != 999999){ // skip fakes, don't store them at all in sim map
      simTkMap[extra.mcTrackID()].push_back(track.label()); // store vector of reco tk labels, mapped to the sim track label (i.e. mcTrackID)
    }
  }

  for (auto&& simTkMatches : simTkMap){
    if (simTkMatches.second.size() < 2) { // no duplicates
      auto& extra(evt_extras[simTkMatches.second[0]]);
      extra.setMCDuplicateInfo(0,bool(false));
    }
    else{ // sort duplicates (ghosts) to keep best one --> most hits, lowest chi2

      // really should sort on indices with a reduced data structure... this is a hack way to do this for now...
      // e.g. std::tuple<int, int, float>, (label, nHits, chi2)
      TrackVec tmpMatches;
      for (auto&& label : simTkMatches.second) { // loop over vector of reco track labels, push back the track with each label 
	tmpMatches.push_back(evt_tracks[label]);
      }
      std::sort(tmpMatches.begin(), tmpMatches.end(), sortByHitsChi2); // sort the tracks
      for (auto itrack = 0; itrack < tmpMatches.size(); itrack++){ // loop over sorted tracks, now make the vector of sorted labels match
 	simTkMatches.second[itrack] = tmpMatches[itrack].label();
      }
      
      int duplicateID = 0;
      for (auto&& label : simTkMatches.second){ // loop over vector of reco tracsk 
        auto& extra(evt_extras[label]);
        extra.setMCDuplicateInfo(duplicateID,bool(true));
        duplicateID++; // used in fake rate trees!
      } 
    }
  }
}

void TTreeValidation::makeSeedTkToRecoTkMaps(Event& ev){
  std::lock_guard<std::mutex> locker(glock_); 

  // map seed to reco tracks --> seed track collection assumed to map to itself, unless we make some cuts
  mapSeedTkToRecoTk(ev.candidateTracks_,ev.candidateTracksExtra_,seedToBuildMap_);
  mapSeedTkToRecoTk(ev.fitTracks_,ev.fitTracksExtra_,seedToFitMap_);
}

void TTreeValidation::mapSeedTkToRecoTk(const TrackVec& evt_tracks, const TrackExtraVec& evt_extras, TkIDToTkIDMap& seedTkMap){
  for (auto&& track : evt_tracks){
    seedTkMap[evt_extras[track.label()].seedID()] = track.label();
  }
}

void TTreeValidation::fillSegmentTree(const BinInfoMap& segmentMap, int evtID){
  std::lock_guard<std::mutex> locker(glock_);

  evtID_seg_ = evtID;
  for (int i = 0; i < Config::nLayers; i++) {
    layer_seg_ = i;
    for (int j = 0; j < Config::nEtaPart; j++) {
      etabin_seg_ = j;
      for (int k = 0; k < Config::nPhiPart; k++) {
	phibin_seg_ = k;
	nHits_seg_  = segmentMap[i][j][k].second;

	segtree_->Fill();
      }
    }
  }
}

void TTreeValidation::fillBranchTree(int evtID)
{
  std::lock_guard<std::mutex> locker(glock_);
  
  evtID_br_  = evtID;
  for (TkIDToBVVMMIter seediter = seedToBranchValVecLayMapMap_.begin(); seediter != seedToBranchValVecLayMapMap_.end(); ++seediter){
    seedID_br_ = (*seediter).first;
    for (BVVLMiter layiter = (*seediter).second.begin(); layiter != (*seediter).second.end(); ++layiter){
      const auto& BranchValVec((*layiter).second);
      const int cands = BranchValVec.size();
      layer_  = (*layiter).first; // first index here is layer

      // clear vectors before filling
      candEtaPhiBins_.clear();
      candHits_.clear();
      candBranches_.clear();
      candnSigmaDeta_.clear();
      candnSigmaDphi_.clear();
        
      // totals
      std::vector<int> candEtaPhiBins(cands);
      std::vector<int> candHits(cands);
      std::vector<int> candBranches(cands);

      // unique hits, etaphibins, branches...
      std::unordered_map<int, bool> uniqueEtaPhiBins; // once a bin, branch, hit is used, set to true to count it only once. take size of map as "uniques"
      std::unordered_map<int, bool> uniqueHits;
      std::unordered_map<int, bool> uniqueBranches;
  
      // nSigmaDeta/phi vec
      std::vector<float> candnSigmaDeta(cands);
      std::vector<float> candnSigmaDphi(cands);

      for (int cand = 0; cand < cands; cand++){ // loop over input candidates at this layer for this seed
	const auto& BranchVal(BranchValVec[cand]); // grab the branch validation object

	////////////////////////////////////
	//  EtaPhiBins Explored Counting  //
	////////////////////////////////////
	
	// set values for total etaphibins explored per candidate, also unique ones for all candidates
	if (BranchVal.phiBinPlus >= BranchVal.phiBinMinus){ // count the number of eta/phi bins explored if no phi wrapping
	  candEtaPhiBins[cand] = (BranchVal.etaBinPlus-BranchVal.etaBinMinus+1)*(BranchVal.phiBinPlus-BranchVal.phiBinMinus+1); // total etaphibins for this candidate

	  // no phi wrap to count uniques
	  for (int ibin = BranchVal.phiBinMinus; ibin <= BranchVal.phiBinPlus; ibin++){
	    uniqueEtaPhiBins[ibin] = true;
	  }
	}
	else{ // if phi wrapping, count need to do this obnoxious counting
	  candEtaPhiBins[cand] = (BranchVal.etaBinPlus-BranchVal.etaBinMinus+1)*(Config::nPhiPart-BranchVal.phiBinMinus+BranchVal.phiBinPlus+1); // total etaphibins for this candidate with phi wrapping

	  // use phi wrapping to count uniques
	  for (int ibin = BranchVal.phiBinMinus; ibin < Config::nPhiPart; ibin++){
	    uniqueEtaPhiBins[ibin] = true;
	  }
	  for (int ibin = 0; ibin <= BranchVal.phiBinPlus; ibin++){
	    uniqueEtaPhiBins[ibin] = true;
	  }
	}
	
	//////////////////////////////
	//  Hits Explored Counting  //
	//////////////////////////////

	candHits[cand] = BranchVal.cand_hit_indices.size(); // set value of nHits explored per input cand for this seed+layer
	for (auto&& cand_hit_idx : BranchVal.cand_hit_indices){ // save unique hits 
	  uniqueHits[cand_hit_idx] = true;
	}

	/////////////////////////////////
	//  Branches Created Counting  //
	/////////////////////////////////

	candBranches[cand] = BranchVal.branch_hit_indices.size(); // set values for total branches created per input candidate before chi2 max cand cleaning (for this seed+layer)
	for (auto&& branch_hit_idx : BranchVal.branch_hit_indices){ // save unique branches --> only consider adding one ghost temp per input?
	  uniqueBranches[branch_hit_idx] = true;
	}

	// store the parameters for width of eta/phi windows for this input candidates building

	candnSigmaDeta[cand] = BranchVal.nSigmaDeta;
	candnSigmaDphi[cand] = BranchVal.nSigmaDphi;
      } // end loop over input candidates for given seed/layer
    
      // fill the rest of the tree after identifying uniques

      cands_  = cands;
      
      // will use these to create total plots (once summed by looping over the vectors), and make per input candidate by looping over these vectors individually
      
      candEtaPhiBins_ = candEtaPhiBins;
      candHits_       = candHits;
      candBranches_   = candBranches;

      // will use this directly to make unique plots by just taking the size of these unordered maps with "trues"

      uniqueEtaPhiBins_ = uniqueEtaPhiBins.size();
      uniqueHits_       = uniqueHits.size();
      uniqueBranches_   = uniqueBranches.size();
      
      // will normalize to equal area, loop over vectors to plot directly

      candnSigmaDeta_ = candnSigmaDeta;
      candnSigmaDphi_ = candnSigmaDphi;
    
      tree_br_->Fill();  // fill once per layer per seed
    } // end loop over layers
  } // end loop over seeds
}

void TTreeValidation::fillEfficiencyTree(const Event& ev){
  std::lock_guard<std::mutex> locker(glock_);

  auto ievt = ev.evtID();
  auto& evt_sim_tracks   = ev.simTracks_;
  auto& evt_seed_tracks  = ev.seedTracks_;
  auto& evt_seed_extras  = ev.seedTracksExtra_;
  auto& evt_build_tracks = ev.candidateTracks_;
  auto& evt_build_extras = ev.candidateTracksExtra_;
  auto& evt_fit_tracks   = ev.fitTracks_;
  auto& evt_fit_extras   = ev.fitTracksExtra_;

  for (auto&& simtrack : evt_sim_tracks){
    evtID_eff_ = ievt;
    mcID_eff_  = simtrack.label();

    // generated values
    pt_mc_gen_eff_  = simtrack.pT(); 
    phi_mc_gen_eff_ = simtrack.momPhi();
    eta_mc_gen_eff_ = simtrack.momEta();
    nHits_mc_eff_   = simtrack.nFoundHits(); // could be that the sim track skips layers!

    // matched seed track
    if (simToSeedMap_.count(mcID_eff_)){ // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToSeedMap_[matched SimID][first element in vector]
      auto& seedtrack = evt_seed_tracks[simToSeedMap_[mcID_eff_][0]]; // returns seedTrack best matched to sim track
      auto& seedextra = evt_seed_extras[seedtrack.label()]; // returns track extra best aligned with seed track
      mcmask_seed_eff_ = 1; // quick logic for matched

      seedID_seed_eff_ = seedextra.seedID(); 

      // use this to access correct sim track layer params
      const float layer = seedtrack.foundLayers().back(); // last layer seed ended up on
      const TrackState & initLayTS = simTkTSVecMap_[mcID_eff_][layer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps

      pt_mc_seed_eff_  = initLayTS.pT();
      phi_mc_seed_eff_ = initLayTS.momPhi();
      eta_mc_seed_eff_ = initLayTS.momEta();

      pt_seed_eff_   = seedtrack.pT();
      ept_seed_eff_  = seedtrack.epT();
      phi_seed_eff_  = seedtrack.momPhi();
      ephi_seed_eff_ = seedtrack.emomPhi();
      eta_seed_eff_  = seedtrack.momEta();
      eeta_seed_eff_ = seedtrack.emomEta();

      // rest of mc info
      nHits_seed_eff_           = seedtrack.nFoundHits();
      nHitsMatched_seed_eff_    = seedextra.nHitsMatched();
      fracHitsMatched_seed_eff_ = float(nHitsMatched_seed_eff_) / float(nHits_seed_eff_);

      hitchi2_seed_eff_   = -10; //seedtrack.chi2(); // currently not being used
      helixchi2_seed_eff_ = computeHelixChi2(initLayTS.parameters,seedtrack.parameters(),seedtrack.errors());

      duplmask_seed_eff_   = seedextra.isDuplicate(); 
      nTkMatches_seed_eff_ = simToSeedMap_[mcID_eff_].size(); // n reco matches to this sim track.
    }
    else{ // unmatched simTracks ... put -99 for all reco values to denote unmatched
      mcmask_seed_eff_  = 0; // quick logic for not matched
      
      seedID_seed_eff_ = -99;
      
      pt_mc_seed_eff_  = -99;
      phi_mc_seed_eff_ = -99;
      eta_mc_seed_eff_ = -99;

      pt_seed_eff_   = -99;
      ept_seed_eff_  = -99;
      phi_seed_eff_  = -99;
      ephi_seed_eff_ = -99;
      eta_seed_eff_  = -99;
      eeta_seed_eff_ = -99;

      nHits_seed_eff_           = -99;
      nHitsMatched_seed_eff_    = -99;
      fracHitsMatched_seed_eff_ = -99;
 
      hitchi2_seed_eff_   = -99;
      helixchi2_seed_eff_ = -99;
      
      duplmask_seed_eff_   = 2; // mask means unmatched sim track
      nTkMatches_seed_eff_ = -99; // unmatched
    }

    // matched build track
    if (simToBuildMap_.count(mcID_eff_)){ // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToBuildMap_[matched SimID][first element in vector]
      auto& buildtrack = evt_build_tracks[simToBuildMap_[mcID_eff_][0]]; // returns buildTrack best matched to sim track
      auto& buildextra = evt_build_extras[buildtrack.label()]; // returns track extra best aligned with build track
      mcmask_build_eff_ = 1; // quick logic for matched

      seedID_build_eff_ = buildextra.seedID(); 

      // use this to access correct sim track layer params
      const float layer = buildtrack.foundLayers().back(); // last layer build ended up on
      const TrackState & initLayTS = simTkTSVecMap_[mcID_eff_][layer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps

      pt_mc_build_eff_  = initLayTS.pT();
      phi_mc_build_eff_ = initLayTS.momPhi();
      eta_mc_build_eff_ = initLayTS.momEta();

      pt_build_eff_   = buildtrack.pT();
      ept_build_eff_  = buildtrack.epT();
      phi_build_eff_  = buildtrack.momPhi();
      ephi_build_eff_ = buildtrack.emomPhi();
      eta_build_eff_  = buildtrack.momEta();
      eeta_build_eff_ = buildtrack.emomEta();
      
      nHits_build_eff_           = buildtrack.nFoundHits();
      nHitsMatched_build_eff_    = buildextra.nHitsMatched();
      fracHitsMatched_build_eff_ = float(nHitsMatched_build_eff_) / float(nHits_build_eff_);

      hitchi2_build_eff_   = buildtrack.chi2(); 
      helixchi2_build_eff_ = computeHelixChi2(initLayTS.parameters,buildtrack.parameters(),buildtrack.errors());

      duplmask_build_eff_   = buildextra.isDuplicate(); 
      nTkMatches_build_eff_ = simToBuildMap_[mcID_eff_].size(); // n reco matches to this sim track.
    }
    else{ // unmatched simTracks ... put -99 for all reco values to denote unmatched
      mcmask_build_eff_  = 0; // quick logic for not matched

      seedID_build_eff_ = -99;

      pt_mc_build_eff_  = -99;
      phi_mc_build_eff_ = -99;
      eta_mc_build_eff_ = -99;

      pt_build_eff_   = -99;
      ept_build_eff_  = -99;
      phi_build_eff_  = -99;
      ephi_build_eff_ = -99;
      eta_build_eff_  = -99;
      eeta_build_eff_ = -99;

      nHits_build_eff_           = -99;
      nHitsMatched_build_eff_    = -99;
      fracHitsMatched_build_eff_ = -99;
 
      hitchi2_build_eff_   = -99;
      helixchi2_build_eff_ = -99;
      
      duplmask_build_eff_   = 2; // mask means unmatched sim track
      nTkMatches_build_eff_ = -99; // unmatched
    }
    
    // matched fit track
    if (simToFitMap_.count(mcID_eff_)){ // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToFitMap_[matched SimID][first element in vector]
      auto& fittrack = evt_fit_tracks[simToFitMap_[mcID_eff_][0]]; // returns fitTrack best matched to sim track
      auto& fitextra = evt_fit_extras[fittrack.label()]; // returns track extra best aligned with fit track
      mcmask_fit_eff_ = 1; // quick logic for matched

      seedID_fit_eff_ = fitextra.seedID(); 

      // use this to access correct sim track layer params

#ifdef INWARDFIT
      const float layer = fittrack.foundLayers().front(); // last layer fit ended up on
#else
      const float layer = fittrack.foundLayers().back(); // last layer fit ended up on
#endif
      const TrackState & initLayTS = simTkTSVecMap_[mcID_eff_][layer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps

      pt_mc_fit_eff_  = initLayTS.pT();
      phi_mc_fit_eff_ = initLayTS.momPhi();
      eta_mc_fit_eff_ = initLayTS.momEta();

      pt_fit_eff_   = fittrack.pT();
      ept_fit_eff_  = fittrack.epT();
      phi_fit_eff_  = fittrack.momPhi();
      ephi_fit_eff_ = fittrack.emomPhi();
      eta_fit_eff_  = fittrack.momEta();
      eeta_fit_eff_ = fittrack.emomEta();
      
      // rest of mc info
      nHits_fit_eff_           = fittrack.nFoundHits();
      nHitsMatched_fit_eff_    = fitextra.nHitsMatched();
      fracHitsMatched_fit_eff_ = float(nHitsMatched_fit_eff_) / float(nHits_fit_eff_);

      hitchi2_fit_eff_   = -10; //fittrack.chi2(); // currently not being used
      helixchi2_fit_eff_ = computeHelixChi2(initLayTS.parameters,fittrack.parameters(),fittrack.errors());

      duplmask_fit_eff_   = fitextra.isDuplicate(); 
      nTkMatches_fit_eff_ = simToFitMap_[mcID_eff_].size(); // n reco matches to this sim track.
    }
    else{ // unmatched simTracks ... put -99 for all reco values to denote unmatched
      mcmask_fit_eff_  = 0; // quick logic for not matched

      seedID_fit_eff_ = -99;

      pt_mc_fit_eff_  = -99;
      phi_mc_fit_eff_ = -99;
      eta_mc_fit_eff_ = -99;

      pt_fit_eff_   = -99;
      ept_fit_eff_  = -99;
      phi_fit_eff_  = -99;
      ephi_fit_eff_ = -99;
      eta_fit_eff_  = -99;
      eeta_fit_eff_ = -99;

      nHits_fit_eff_           = -99;
      nHitsMatched_fit_eff_    = -99;
      fracHitsMatched_fit_eff_ = -99;
 
      hitchi2_fit_eff_   = -99;
      helixchi2_fit_eff_ = -99;
      
      duplmask_fit_eff_   = 2; // mask means unmatched sim track
      nTkMatches_fit_eff_ = -99; // unmatched
    }

    efftree_->Fill(); // fill it once per sim track!
  }
}

void TTreeValidation::fillFakeRateTree(const Event& ev){
  std::lock_guard<std::mutex> locker(glock_);

  auto ievt = ev.evtID();
  auto& evt_seed_tracks  = ev.seedTracks_;
  auto& evt_seed_extras  = ev.seedTracksExtra_;
  auto& evt_build_tracks = ev.candidateTracks_;
  auto& evt_build_extras = ev.candidateTracksExtra_;
  auto& evt_fit_tracks   = ev.fitTracks_;
  auto& evt_fit_extras   = ev.fitTracksExtra_;
  
  for (auto&& seedtrack : evt_seed_tracks){
    evtID_FR_       = ievt;
    auto& seedextra = evt_seed_extras[seedtrack.label()];
    seedID_FR_      = seedextra.seedID();

    // seed info
    seedmask_seed_FR_ = 1; // automatically set to 1, because at the moment no cuts on seeds after conformal+KF fit.  seed triplets filtered by RZ chi2 before fitting. 

    pt_seed_FR_   = seedtrack.pT();
    ept_seed_FR_  = seedtrack.epT();
    phi_seed_FR_  = seedtrack.momPhi();
    ephi_seed_FR_ = seedtrack.emomPhi();
    eta_seed_FR_  = seedtrack.momEta();
    eeta_seed_FR_ = seedtrack.emomEta();

    nHits_seed_FR_           = seedtrack.nFoundHits();
    nHitsMatched_seed_FR_    = seedextra.nHitsMatched();
    fracHitsMatched_seed_FR_ = float(nHitsMatched_seed_FR_) / float(nHits_seed_FR_);

    hitchi2_seed_FR_ = -10; // seedtrack.chi2(); --> not currently used

    // sim info for seed track
    mcID_seed_FR_ = seedextra.mcTrackID();
    if (mcID_seed_FR_ != 999999){ // store sim info at that final layer!!! --> gen info stored only in eff tree
      mcmask_seed_FR_ = 1; // matched track to sim

      const float layer = seedtrack.foundLayers().back(); // last layer fit ended up on
      const TrackState & initLayTS = simTkTSVecMap_[mcID_seed_FR_][layer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps

      pt_mc_seed_FR_    = initLayTS.pT();
      phi_mc_seed_FR_   = initLayTS.momPhi();
      eta_mc_seed_FR_   = initLayTS.momEta();
      nHits_mc_seed_FR_ = simTkTSVecMap_[mcID_seed_FR_].size();

      helixchi2_seed_FR_ = computeHelixChi2(initLayTS.parameters,seedtrack.parameters(),seedtrack.errors());

      duplmask_seed_FR_   = seedextra.isDuplicate();
      iTkMatches_seed_FR_ = seedextra.duplicateID(); // ith duplicate seed track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"      
    }
    else{
      mcmask_seed_FR_ = 0;   // fake track (unmatched track)
          
      // -99 for all sim info for reco tracks not associated to reco tracks
      pt_mc_seed_FR_  = -99;
      phi_mc_seed_FR_ = -99;
      eta_mc_seed_FR_ = -99;
      nHits_mc_seed_FR_ = -99;

      duplmask_seed_FR_   = 2; // see notation above      
      iTkMatches_seed_FR_ = -99;  
    }

    //==========================//
    
    // fill build information if track still alive
    if (seedToBuildMap_.count(seedID_FR_)){
      seedmask_build_FR_ = 1; // quick logic

      auto& buildtrack = evt_build_tracks[seedToBuildMap_[seedID_FR_]];
      auto& buildextra = evt_build_extras[buildtrack.label()];

      pt_build_FR_   = buildtrack.pT();
      ept_build_FR_  = buildtrack.epT();
      phi_build_FR_  = buildtrack.momPhi();
      ephi_build_FR_ = buildtrack.emomPhi();
      eta_build_FR_  = buildtrack.momEta();
      eeta_build_FR_ = buildtrack.emomEta();

      nHits_build_FR_           = buildtrack.nFoundHits();
      nHitsMatched_build_FR_    = buildextra.nHitsMatched();
      fracHitsMatched_build_FR_ = float(nHitsMatched_build_FR_) / float(nHits_build_FR_);

      hitchi2_build_FR_ = buildtrack.chi2();

      // sim info for build track
      mcID_build_FR_  = buildextra.mcTrackID();
      if (mcID_build_FR_ != 999999){ // build track matched to seed and sim 
	mcmask_build_FR_ = 1; // matched track to sim
	
	const float layer = buildtrack.foundLayers().back(); // last layer fit ended up on
	const TrackState & initLayTS = simTkTSVecMap_[mcID_build_FR_][layer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps
	
	pt_mc_build_FR_    = initLayTS.pT();
	phi_mc_build_FR_   = initLayTS.momPhi();
	eta_mc_build_FR_   = initLayTS.momEta();
	nHits_mc_build_FR_ = simTkTSVecMap_[mcID_build_FR_].size();
	
	helixchi2_build_FR_ = computeHelixChi2(initLayTS.parameters,buildtrack.parameters(),buildtrack.errors());

	duplmask_build_FR_   = buildextra.isDuplicate();
	iTkMatches_build_FR_ = buildextra.duplicateID(); // ith duplicate build track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"      
      }
      else{ // build track matched only to seed not to sim
	mcmask_build_FR_ = 0;   // fake track (unmatched track)
	
	// -99 for all sim info for reco tracks not associated to reco tracks
	pt_mc_build_FR_  = -99;
	phi_mc_build_FR_ = -99;
	eta_mc_build_FR_ = -99;
	nHits_mc_build_FR_ = -99;

	duplmask_build_FR_   = 2;
	iTkMatches_build_FR_ = -99;
      } // matched seed to build, not build to sim
    }

    else { // seed has no matching build track (therefore no matching sim to build track)
      seedmask_build_FR_ = 0; // quick logic

      // -100 for all reco info as no actual build track for this seed
      pt_build_FR_   = -100;
      ept_build_FR_  = -100;
      phi_build_FR_  = -100;
      ephi_build_FR_ = -100;
      eta_build_FR_  = -100;
      eeta_build_FR_ = -100;
      
      nHits_build_FR_ = -100;
      nHitsMatched_build_FR_ = -100;
      fracHitsMatched_build_FR_ = -100;

      hitchi2_build_FR_  = -100; 
      
      // keep -100 for all sim variables as no such reco exists for this seed
      mcmask_build_FR_ = -1; // do not want to count towards build FR
      mcID_build_FR_   = -100;
	
      pt_mc_build_FR_  = -100;
      phi_mc_build_FR_ = -100;
      eta_mc_build_FR_ = -100;
      nHits_mc_build_FR_ = -100;

      helixchi2_build_FR_ = -100;

      duplmask_build_FR_   = -1;
      iTkMatches_build_FR_ = -100;
    }

    //============================// fit tracks
    if (seedToFitMap_.count(seedID_FR_)){
      seedmask_fit_FR_ = 1; // quick logic

      auto& fittrack = evt_fit_tracks[seedToFitMap_[seedID_FR_]];
      auto& fitextra = evt_fit_extras[fittrack.label()];

      pt_fit_FR_   = fittrack.pT();
      ept_fit_FR_  = fittrack.epT();
      phi_fit_FR_  = fittrack.momPhi();
      ephi_fit_FR_ = fittrack.emomPhi();
      eta_fit_FR_  = fittrack.momEta();
      eeta_fit_FR_ = fittrack.emomEta();

      nHits_fit_FR_           = fittrack.nFoundHits();
      nHitsMatched_fit_FR_    = fitextra.nHitsMatched();
      fracHitsMatched_fit_FR_ = float(nHitsMatched_fit_FR_) / float(nHits_fit_FR_);

      hitchi2_fit_FR_ = -10; //fittrack.chi2() --> currently not used

      // sim info for fit track
      mcID_fit_FR_  = fitextra.mcTrackID();
      if (mcID_fit_FR_ != 999999){ // fit track matched to seed and sim 
	mcmask_fit_FR_ = 1; // matched track to sim

#ifdef INWARDFIT
	const float layer = fittrack.foundLayers().front(); // last layer fiting ended up on
#else
	const float layer = fittrack.foundLayers().back(); // last layer fiting ended up on
#endif
	const TrackState & initLayTS = simTkTSVecMap_[mcID_fit_FR_][layer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps
   
	pt_mc_fit_FR_    = initLayTS.pT(); // again, pt of mc truth at the layer the seed ends
	phi_mc_fit_FR_   = initLayTS.momPhi();
	eta_mc_fit_FR_   = initLayTS.momEta();
	nHits_mc_fit_FR_ = simTkTSVecMap_[mcID_fit_FR_].size();

	helixchi2_fit_FR_ = computeHelixChi2(initLayTS.parameters,fittrack.parameters(),fittrack.errors());

	duplmask_fit_FR_   = fitextra.isDuplicate();
	iTkMatches_fit_FR_ = fitextra.duplicateID(); // ith duplicate fit track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"
      }
      else{ // fit track matched only to seed not to sim
	mcmask_fit_FR_ = 0;   // fake track (unmatched track)
	
	// -99 for all sim info for reco tracks not associated to reco tracks
	pt_mc_fit_FR_  = -99;
	phi_mc_fit_FR_ = -99;
	eta_mc_fit_FR_ = -99;
	nHits_mc_fit_FR_ = -99;

	duplmask_fit_FR_   = 2;
	iTkMatches_fit_FR_ = -99;
      } // matched seed to fit, not fit to sim
    }

    else { // seed has no matching fit track (therefore no matching sim to fit track)
      seedmask_fit_FR_ = 0; // quick logic

      // -100 for all reco info as no actual fit track for this seed
      pt_fit_FR_   = -100;
      ept_fit_FR_  = -100;
      phi_fit_FR_  = -100;
      ephi_fit_FR_ = -100;
      eta_fit_FR_  = -100;
      eeta_fit_FR_ = -100;
      
      nHits_fit_FR_ = -100;
      nHitsMatched_fit_FR_ = -100;
      fracHitsMatched_fit_FR_ = -100;

      hitchi2_fit_FR_  = -100; 
      
      // keep -100 for all sim variables as no such reco exists for this seed
      mcmask_fit_FR_ = -1; // do not want to count towards fit FR
      mcID_fit_FR_   = -100;
	
      pt_mc_fit_FR_  = -100;
      phi_mc_fit_FR_ = -100;
      eta_mc_fit_FR_ = -100;
      nHits_mc_fit_FR_ = -100;

      helixchi2_fit_FR_ = -100;

      duplmask_fit_FR_   = -1;
      iTkMatches_fit_FR_ = -100;
    }
        
    fakeratetree_->Fill(); // fill once per seed!
  }// end of seed to seed loop
}

void TTreeValidation::fillGeometryTree(const Event& ev){
  std::lock_guard<std::mutex> locker(glock_);

  auto ievt = ev.evtID();
  auto& evt_sim_tracks   = ev.simTracks_;
  auto& evt_seed_tracks  = ev.seedTracks_;
  auto& evt_seed_extras  = ev.seedTracksExtra_;
  auto& evt_fit_tracks   = ev.fitTracks_;
  auto& evt_fit_extras   = ev.fitTracksExtra_;
  auto& evt_lay_hits     = ev.layerHits_;

  for (auto&& simtrack : evt_sim_tracks){
    evtID_geo_ = ievt;
    mcID_geo_  = simtrack.label();

    // for Beamspot plots
    x_mc_gen_vrx_geo_ = simtrack.x();
    y_mc_gen_vrx_geo_ = simtrack.y();
    z_mc_gen_vrx_geo_ = simtrack.z();

    // clear for detector sim plots
    x_mc_reco_hit_geo_.clear();
    y_mc_reco_hit_geo_.clear();
    z_mc_reco_hit_geo_.clear();
    for (auto&& simhit : simtrack.hitsVector(evt_lay_hits)){ // assume one hit per layer
      x_mc_reco_hit_geo_.push_back(simhit.x());
      y_mc_reco_hit_geo_.push_back(simhit.y());
      z_mc_reco_hit_geo_.push_back(simhit.z());
    }

    // clear vectors for position pulls
    layers_seed_geo_.clear();
    x_lay_seed_geo_.clear();
    y_lay_seed_geo_.clear();
    z_lay_seed_geo_.clear();
    ex_lay_seed_geo_.clear();
    ey_lay_seed_geo_.clear();
    ez_lay_seed_geo_.clear();

    layers_fit_geo_.clear();
    x_lay_fit_geo_.clear();
    y_lay_fit_geo_.clear();
    z_lay_fit_geo_.clear();
    ex_lay_fit_geo_.clear();
    ey_lay_fit_geo_.clear();
    ez_lay_fit_geo_.clear();

    // matched seed track
    if (simToSeedMap_.count(mcID_geo_)){ // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToSeedMap_[matched SimID][first element in vector]
      auto& seedtrack = evt_seed_tracks[simToSeedMap_[mcID_geo_][0]]; // returns seedTrack best matched to sim track
      auto& seedextra = evt_seed_extras[seedtrack.label()]; // returns track extra best aligned with seed track
      mcmask_seed_geo_ = 1; // quick logic for matched

      seedID_seed_geo_ = seedextra.seedID(); 

      // position pull info
      const TSLayerPairVec & seedTSLayerPairVec = seedTkTSLayerPairVecMap_[seedID_seed_geo_];
      for (int ilay = 0; ilay < seedTSLayerPairVec.size(); ilay++){ // loop over layers present in pair vector
	layers_seed_geo_.push_back(seedTSLayerPairVec[ilay].first); // want to push back the ACTUAL layers stored, not the index of the loop!
	
	x_lay_seed_geo_.push_back(seedTSLayerPairVec[ilay].second.x()); // therefore, trackstate filled in sync with the layer it was saved on for the vector
	y_lay_seed_geo_.push_back(seedTSLayerPairVec[ilay].second.y());
	z_lay_seed_geo_.push_back(seedTSLayerPairVec[ilay].second.z());
	
	ex_lay_seed_geo_.push_back(seedTSLayerPairVec[ilay].second.exx());
	ey_lay_seed_geo_.push_back(seedTSLayerPairVec[ilay].second.eyy());
	ez_lay_seed_geo_.push_back(seedTSLayerPairVec[ilay].second.ezz());
      }
    }
    else{ // unmatched sim track for any seed track
      layers_seed_geo_.push_back(-1); // mask per layer
      
      x_lay_seed_geo_.push_back(-2000);
      y_lay_seed_geo_.push_back(-2000);
      z_lay_seed_geo_.push_back(-2000);
      
      ex_lay_seed_geo_.push_back(-2000);
      ey_lay_seed_geo_.push_back(-2000);
      ez_lay_seed_geo_.push_back(-2000);
    }

    if (simToFitMap_.count(mcID_geo_)){ // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToFitMap_[matched SimID][first element in vector]
      auto& fittrack = evt_fit_tracks[simToFitMap_[mcID_geo_][0]]; // returns fitTrack best matched to sim track
      auto& fitextra = evt_fit_extras[fittrack.label()]; // returns track extra best aligned with fit track
      mcmask_fit_geo_ = 1; // quick logic for matched

      seedID_fit_geo_ = fitextra.seedID(); 

      // position pull info
      const TSLayerPairVec & fitTSLayerPairVec = fitTkTSLayerPairVecMap_[seedID_fit_geo_];
      for (int ilay = 0; ilay < fitTSLayerPairVec.size(); ilay++){ // loop over layers present in pair vector
	layers_fit_geo_.push_back(fitTSLayerPairVec[ilay].first); // want to push back the ACTUAL layers stored, not the index of the loop!
	
	x_lay_fit_geo_.push_back(fitTSLayerPairVec[ilay].second.x()); // therefore, trackstate filled in sync with the layer it was saved on for the vector
	y_lay_fit_geo_.push_back(fitTSLayerPairVec[ilay].second.y());
	z_lay_fit_geo_.push_back(fitTSLayerPairVec[ilay].second.z());
	
	ex_lay_fit_geo_.push_back(fitTSLayerPairVec[ilay].second.exx());
	ey_lay_fit_geo_.push_back(fitTSLayerPairVec[ilay].second.eyy());
	ez_lay_fit_geo_.push_back(fitTSLayerPairVec[ilay].second.ezz());
      }
    }
    else{
      layers_fit_geo_.push_back(-1); // mask per layer
      
      x_lay_fit_geo_.push_back(-2000);
      y_lay_fit_geo_.push_back(-2000);
      z_lay_fit_geo_.push_back(-2000);
      
      ex_lay_fit_geo_.push_back(-2000);
      ey_lay_fit_geo_.push_back(-2000);
      ez_lay_fit_geo_.push_back(-2000);
    }
    geotree_->Fill(); // fill once per sim track
  }
}

void TTreeValidation::fillConformalTree(const Event& ev){
  std::lock_guard<std::mutex> locker(glock_);

  auto ievt = ev.evtID();
  auto& evt_sim_tracks   = ev.simTracks_;
  auto& evt_seed_tracks  = ev.seedTracks_;
  auto& evt_seed_extras  = ev.seedTracksExtra_;
  auto& evt_fit_tracks   = ev.fitTracks_;
  auto& evt_fit_extras   = ev.fitTracksExtra_;
  auto& evt_lay_hits     = ev.layerHits_;

  for (auto&& simtrack : evt_sim_tracks){
    evtID_cf_ = ievt;
    mcID_cf_  = simtrack.label();

    if (simToSeedMap_.count(mcID_cf_)){ // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToSeedMap_[matched SimID][first element in vector]
      auto& seedtrack = evt_seed_tracks[simToSeedMap_[mcID_cf_][0]]; // returns seedTrack best matched to sim track
      auto& seedextra = evt_seed_extras[seedtrack.label()]; // returns track extra best aligned with seed track
      mcmask_seed_cf_ = 1; // quick logic for matched

      seedID_seed_cf_ = seedextra.seedID(); 

      // Conformal fit stuff to match how it was done before
      const float cflayer = seedtrack.foundLayers().front(); //  layer for which cf parameters are calculated with respect to (either first or last layer of track!)
      const TrackState & cfInitLayTS = simTkTSVecMap_[mcID_cf_][cflayer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps
      x_mc_seed_cf_ = cfInitLayTS.x();
      y_mc_seed_cf_ = cfInitLayTS.y();
      z_mc_seed_cf_ = cfInitLayTS.z();

      px_mc_seed_cf_ = cfInitLayTS.px();
      py_mc_seed_cf_ = cfInitLayTS.py();
      pz_mc_seed_cf_ = cfInitLayTS.pz();

      pt_mc_seed_cf_    = cfInitLayTS.pT();
      invpt_mc_seed_cf_ = cfInitLayTS.invpT();
      phi_mc_seed_cf_   = cfInitLayTS.momPhi();
      theta_mc_seed_cf_ = cfInitLayTS.theta();

      const TrackState & cfSeedTS = seedTkCFMap_[seedID_seed_cf_];

      x_seed_cf_ = cfSeedTS.x();
      y_seed_cf_ = cfSeedTS.y();
      z_seed_cf_ = cfSeedTS.z();
      
      ex_seed_cf_ = cfSeedTS.exx();
      ey_seed_cf_ = cfSeedTS.eyy();
      ez_seed_cf_ = cfSeedTS.ezz();

      px_seed_cf_ = cfSeedTS.px();
      py_seed_cf_ = cfSeedTS.py();
      pz_seed_cf_ = cfSeedTS.pz();

      epx_seed_cf_ = cfSeedTS.epxpx();
      epy_seed_cf_ = cfSeedTS.epypy();
      epz_seed_cf_ = cfSeedTS.epzpz();

      pt_seed_cf_    = cfSeedTS.pT();
      invpt_seed_cf_ = cfSeedTS.invpT();
      phi_seed_cf_   = cfSeedTS.momPhi();
      theta_seed_cf_ = cfSeedTS.theta();

      ept_seed_cf_    = cfSeedTS.epT();
      einvpt_seed_cf_ = cfSeedTS.einvpT();
      ephi_seed_cf_   = cfSeedTS.emomPhi();
      etheta_seed_cf_ = cfSeedTS.etheta();
    }
    else{ // no matched seed form sim
      x_mc_seed_cf_ = -2000;
      y_mc_seed_cf_ = -2000;
      z_mc_seed_cf_ = -2000;

      px_mc_seed_cf_ = -99;
      py_mc_seed_cf_ = -99;
      pz_mc_seed_cf_ = -99;

      pt_mc_seed_cf_    = -99;
      invpt_mc_seed_cf_ = -99;
      phi_mc_seed_cf_   = -99;
      theta_mc_seed_cf_ = -99;

      x_seed_cf_ = -2000;
      y_seed_cf_ = -2000;
      z_seed_cf_ = -2000;
      
      ex_seed_cf_ = -2000;
      ey_seed_cf_ = -2000;
      ez_seed_cf_ = -2000;

      px_seed_cf_ = -99;
      py_seed_cf_ = -99;
      pz_seed_cf_ = -99;

      epx_seed_cf_ = -99;
      epy_seed_cf_ = -99;
      epz_seed_cf_ = -99;

      pt_seed_cf_    = -99;
      invpt_seed_cf_ = -99;
      phi_seed_cf_   = -99;
      theta_seed_cf_ = -99;

      ept_seed_cf_    = -99;
      einvpt_seed_cf_ = -99;
      ephi_seed_cf_   = -99;
      etheta_seed_cf_ = -99;
    }

    if (simToFitMap_.count(mcID_cf_)){ // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToFitMap_[matched SimID][first element in vector]
      auto& fittrack = evt_fit_tracks[simToFitMap_[mcID_cf_][0]]; // returns fitTrack best matched to sim track
      auto& fitextra = evt_fit_extras[fittrack.label()]; // returns track extra best aligned with fit track
      mcmask_fit_cf_ = 1; // quick logic for matched

      seedID_fit_cf_ = fitextra.seedID(); 

      // Conformal fit stuff to match how it was done before
#ifdef INWARDFIT
      const float cflayer = fittrack.foundLayers().back(); // last layer fit ended up on
#else
      const float cflayer = fittrack.foundLayers().front(); // last layer fit ended up on
#endif
      const TrackState & cfInitLayTS = simTkTSVecMap_[mcID_cf_][cflayer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps
      x_mc_fit_cf_ = cfInitLayTS.x();
      y_mc_fit_cf_ = cfInitLayTS.y();
      z_mc_fit_cf_ = cfInitLayTS.z();

      px_mc_fit_cf_ = cfInitLayTS.px();
      py_mc_fit_cf_ = cfInitLayTS.py();
      pz_mc_fit_cf_ = cfInitLayTS.pz();

      pt_mc_fit_cf_    = cfInitLayTS.pT();
      invpt_mc_fit_cf_ = cfInitLayTS.invpT();
      phi_mc_fit_cf_   = cfInitLayTS.momPhi();
      theta_mc_fit_cf_ = cfInitLayTS.theta();

      const TrackState & cfFitTS = fitTkCFMap_[seedID_fit_cf_];

      x_fit_cf_ = cfFitTS.x();
      y_fit_cf_ = cfFitTS.y();
      z_fit_cf_ = cfFitTS.z();
      
      ex_fit_cf_ = cfFitTS.exx();
      ey_fit_cf_ = cfFitTS.eyy();
      ez_fit_cf_ = cfFitTS.ezz();

      px_fit_cf_ = cfFitTS.px();
      py_fit_cf_ = cfFitTS.py();
      pz_fit_cf_ = cfFitTS.pz();

      epx_fit_cf_ = cfFitTS.epxpx();
      epy_fit_cf_ = cfFitTS.epypy();
      epz_fit_cf_ = cfFitTS.epzpz();

      pt_fit_cf_    = cfFitTS.pT();
      invpt_fit_cf_ = cfFitTS.invpT();
      phi_fit_cf_   = cfFitTS.momPhi();
      theta_fit_cf_ = cfFitTS.theta();

      ept_fit_cf_    = cfFitTS.epT();
      einvpt_fit_cf_ = cfFitTS.einvpT();
      ephi_fit_cf_   = cfFitTS.emomPhi();
      etheta_fit_cf_ = cfFitTS.etheta();
    }
    else{ // no matched fit form sim
      x_mc_fit_cf_ = -2000;
      y_mc_fit_cf_ = -2000;
      z_mc_fit_cf_ = -2000;

      px_mc_fit_cf_ = -99;
      py_mc_fit_cf_ = -99;
      pz_mc_fit_cf_ = -99;

      pt_mc_fit_cf_    = -99;
      invpt_mc_fit_cf_ = -99;
      phi_mc_fit_cf_   = -99;
      theta_mc_fit_cf_ = -99;

      x_fit_cf_ = -2000;
      y_fit_cf_ = -2000;
      z_fit_cf_ = -2000;
      
      ex_fit_cf_ = -2000;
      ey_fit_cf_ = -2000;
      ez_fit_cf_ = -2000;

      px_fit_cf_ = -99;
      py_fit_cf_ = -99;
      pz_fit_cf_ = -99;

      epx_fit_cf_ = -99;
      epy_fit_cf_ = -99;
      epz_fit_cf_ = -99;

      pt_fit_cf_    = -99;
      invpt_fit_cf_ = -99;
      phi_fit_cf_   = -99;
      theta_fit_cf_ = -99;

      ept_fit_cf_    = -99;
      einvpt_fit_cf_ = -99;
      ephi_fit_cf_   = -99;
      etheta_fit_cf_ = -99;
    }
    cftree_->Fill();
  }
}

void TTreeValidation::fillConfigTree(const std::vector<double>& ticks){
  simtime_   = ticks[0];
  segtime_   = ticks[1];
  seedtime_  = ticks[2];
  buildtime_ = ticks[3];
  fittime_   = ticks[4];
  hlvtime_   = ticks[5];

  Ntracks_ = Config::nTracks;
  Nevents_ = Config::nEvents;

  nLayers_ = Config::nLayers;
  fRadialSpacing_ = Config::fRadialSpacing;
  fRadialExtent_  = Config::fRadialExtent;
  fInnerSensorSize_ = Config::fInnerSensorSize;
  fOuterSensorSize_ = Config::fOuterSensorSize;
  fEtaDet_ = Config::fEtaDet;

  nPhiPart_   = Config::nPhiPart;
  fPhiFactor_ = Config::fPhiFactor;
  nEtaPart_   = Config::nEtaPart;

  nlayers_per_seed_ = Config::nlayers_per_seed;
  maxCand_ = Config::maxCandsPerSeed;
  chi2Cut_ = Config::chi2Cut;
  nSigma_  = Config::nSigma;
  minDPhi_ = Config::minDPhi;
  maxDPhi_ = Config::maxDPhi;
  minDEta_ = Config::minDEta;
  maxDEta_ = Config::maxDEta;
  
  beamspotX_ = Config::beamspotX;
  beamspotY_ = Config::beamspotY;
  beamspotZ_ = Config::beamspotZ;

  minSimPt_ = Config::minSimPt;
  maxSimPt_ = Config::maxSimPt;

  hitposerrXY_ = Config::hitposerrXY;
  hitposerrZ_  = Config::hitposerrZ;
  hitposerrR_  = Config::hitposerrR;
  varXY_ = Config::varXY;
  varZ_  = Config::varZ;
  
  nTotHit_ = Config::nTotHit;

  ptinverr049_ = Config::ptinverr049;
  phierr049_   = Config::phierr049;
  thetaerr049_ = Config::thetaerr049;
  ptinverr012_ = Config::ptinverr012;
  phierr012_   = Config::phierr012;
  thetaerr012_ = Config::thetaerr012;

  configtree_->Fill();
}

void TTreeValidation::saveTTrees() {  
  std::lock_guard<std::mutex> locker(glock_); 
  f_->cd();
  f_->Write();
  f_->Close();
}             

#endif
