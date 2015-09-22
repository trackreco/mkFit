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
#include "Config.h"
#include "Propagation.h"
#ifndef NO_ROOT

static bool sortByHitsChi2(const Track* cand1, const Track* cand2)
{
  if (cand1->nHits()==cand2->nHits()) return cand1->chi2()<cand2->chi2();
  return cand1->nHits()>cand2->nHits();
}

inline float computeHelixChi2(const SVector6& simParams, const SVector6& recoParams, const SMatrixSym66& recoErrs)
{ 
  float chi2 = 0;
  for(auto i = 0U; i < 6; i++){
    float delta = simParams.At(i) - recoParams.At(i);
    chi2 += (delta*delta) / recoErrs.At(i,i);
  }
  return chi2 / 5;
}

TTreeValidation::TTreeValidation(std::string fileName)
{
  std::lock_guard<std::mutex> locker(glock_);
  gROOT->ProcessLine("#include <vector>");
  f_ = TFile::Open(fileName.c_str(), "recreate");

  // include ALL config ++ real seeding parameters ...
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

  // segment validation
  segtree_ = new TTree("segtree","segtree");
  segtree_->Branch("evtID",&evtID_seg_);
  segtree_->Branch("layer",&layer_seg_);
  segtree_->Branch("etabin",&etabin_seg_);
  segtree_->Branch("phibin",&phibin_seg_);
  segtree_->Branch("nHits",&nHits_seg_);

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
  
  // efficiency validation
  efftree_ = new TTree("efftree","efftree");
  efftree_->Branch("evtID",&evtID_eff_);
  efftree_->Branch("mcID",&mcID_eff_);

  efftree_->Branch("seedID_seed",&seedID_seed_eff_);
  efftree_->Branch("seedID_build",&seedID_build_eff_);
  efftree_->Branch("seedID_fit",&seedID_fit_eff_);

  efftree_->Branch("pt_mc_gen",&pt_mc_gen_eff_);
  efftree_->Branch("phi_mc_gen",&phi_mc_gen_eff_);
  efftree_->Branch("eta_mc_gen",&eta_mc_gen_eff_);
  efftree_->Branch("nHits_mc",&nHits_mc_eff_);
  
  efftree_->Branch("x_mc_reco_hit",&x_mc_reco_hit_eff_);
  efftree_->Branch("y_mc_reco_hit",&y_mc_reco_hit_eff_);
  efftree_->Branch("z_mc_reco_hit",&z_mc_reco_hit_eff_);
  
  efftree_->Branch("x_mc_gen_vrx",&x_mc_gen_vrx_eff_);
  efftree_->Branch("y_mc_gen_vrx",&y_mc_gen_vrx_eff_);
  efftree_->Branch("z_mc_gen_vrx",&z_mc_gen_vrx_eff_);

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

  // eventually move all conformal stuff to a dedicated separate "conftree" --> separate from efficiency tree as too large

  efftree_->Branch("x_mc_cf_seed",&x_mc_cf_seed_eff_);
  efftree_->Branch("x_cf_seed",&x_cf_seed_eff_);
  efftree_->Branch("ex_cf_seed",&ex_cf_seed_eff_);
  efftree_->Branch("x_mc_cf_fit",&x_mc_cf_fit_eff_);
  efftree_->Branch("x_cf_fit",&x_cf_fit_eff_);
  efftree_->Branch("ex_cf_fit",&ex_cf_fit_eff_);

  efftree_->Branch("y_mc_cf_seed",&y_mc_cf_seed_eff_);
  efftree_->Branch("y_cf_seed",&y_cf_seed_eff_);
  efftree_->Branch("ey_cf_seed",&ey_cf_seed_eff_);
  efftree_->Branch("y_mc_cf_fit",&y_mc_cf_fit_eff_);
  efftree_->Branch("y_cf_fit",&y_cf_fit_eff_);
  efftree_->Branch("ey_cf_fit",&ey_cf_fit_eff_);

  efftree_->Branch("z_mc_cf_seed",&z_mc_cf_seed_eff_);
  efftree_->Branch("z_cf_seed",&z_cf_seed_eff_);
  efftree_->Branch("ez_cf_seed",&ez_cf_seed_eff_);
  efftree_->Branch("z_mc_cf_fit",&z_mc_cf_fit_eff_);
  efftree_->Branch("z_cf_fit",&z_cf_fit_eff_);
  efftree_->Branch("ez_cf_fit",&ez_cf_fit_eff_);

  efftree_->Branch("px_mc_cf_seed",&px_mc_cf_seed_eff_);
  efftree_->Branch("px_cf_seed",&px_cf_seed_eff_);
  efftree_->Branch("epx_cf_seed",&epx_cf_seed_eff_);
  efftree_->Branch("px_mc_cf_fit",&px_mc_cf_fit_eff_);
  efftree_->Branch("px_cf_fit",&px_cf_fit_eff_);
  efftree_->Branch("epx_cf_fit",&epx_cf_fit_eff_);

  efftree_->Branch("py_mc_cf_seed",&py_mc_cf_seed_eff_);
  efftree_->Branch("py_cf_seed",&py_cf_seed_eff_);
  efftree_->Branch("epy_cf_seed",&epy_cf_seed_eff_);
  efftree_->Branch("py_mc_cf_fit",&py_mc_cf_fit_eff_);
  efftree_->Branch("py_cf_fit",&py_cf_fit_eff_);
  efftree_->Branch("epy_cf_fit",&epy_cf_fit_eff_);

  efftree_->Branch("pz_mc_cf_seed",&pz_mc_cf_seed_eff_);
  efftree_->Branch("pz_cf_seed",&pz_cf_seed_eff_);
  efftree_->Branch("epz_cf_seed",&epz_cf_seed_eff_);
  efftree_->Branch("pz_mc_cf_fit",&pz_mc_cf_fit_eff_);
  efftree_->Branch("pz_cf_fit",&pz_cf_fit_eff_);
  efftree_->Branch("epz_cf_fit",&epz_cf_fit_eff_);

  efftree_->Branch("pt_mc_cf_seed",&pt_mc_cf_seed_eff_);
  efftree_->Branch("pt_cf_seed",&pt_cf_seed_eff_);
  efftree_->Branch("ept_cf_seed",&ept_cf_seed_eff_);
  efftree_->Branch("pt_mc_cf_fit",&pt_mc_cf_fit_eff_);
  efftree_->Branch("pt_cf_fit",&pt_cf_fit_eff_);
  efftree_->Branch("ept_cf_fit",&ept_cf_fit_eff_);

  efftree_->Branch("invpt_mc_cf_seed",&invpt_mc_cf_seed_eff_);
  efftree_->Branch("invpt_cf_seed",&invpt_cf_seed_eff_);
  efftree_->Branch("einvpt_cf_seed",&einvpt_cf_seed_eff_);
  efftree_->Branch("invpt_mc_cf_fit",&invpt_mc_cf_fit_eff_);
  efftree_->Branch("invpt_cf_fit",&invpt_cf_fit_eff_);
  efftree_->Branch("einvpt_cf_fit",&einvpt_cf_fit_eff_);

  efftree_->Branch("phi_mc_cf_seed",&phi_mc_cf_seed_eff_);
  efftree_->Branch("phi_cf_seed",&phi_cf_seed_eff_);
  efftree_->Branch("ephi_cf_seed",&ephi_cf_seed_eff_);
  efftree_->Branch("phi_mc_cf_fit",&phi_mc_cf_fit_eff_);
  efftree_->Branch("phi_cf_fit",&phi_cf_fit_eff_);
  efftree_->Branch("ephi_cf_fit",&ephi_cf_fit_eff_);
  
  efftree_->Branch("theta_mc_cf_seed",&theta_mc_cf_seed_eff_);
  efftree_->Branch("theta_cf_seed",&theta_cf_seed_eff_);
  efftree_->Branch("etheta_cf_seed",&etheta_cf_seed_eff_);
  efftree_->Branch("theta_mc_cf_fit",&theta_mc_cf_fit_eff_);
  efftree_->Branch("theta_cf_fit",&theta_cf_fit_eff_);
  efftree_->Branch("etheta_cf_fit",&etheta_cf_fit_eff_);

  // end of conformal stuff

  // position pull info --> should be on it's own tree or not included at all...
  efftree_->Branch("layers_seed",&layers_seed_eff_);
  efftree_->Branch("layers_fit",&layers_fit_eff_);

  efftree_->Branch("x_lay_seed",&x_lay_seed_eff_);
  efftree_->Branch("y_lay_seed",&y_lay_seed_eff_);
  efftree_->Branch("z_lay_seed",&z_lay_seed_eff_);
  efftree_->Branch("ex_lay_seed",&ex_lay_seed_eff_);
  efftree_->Branch("ey_lay_seed",&ey_lay_seed_eff_);
  efftree_->Branch("ez_lay_seed",&ez_lay_seed_eff_);

  efftree_->Branch("x_lay_fit",&x_lay_fit_eff_);
  efftree_->Branch("y_lay_fit",&y_lay_fit_eff_);
  efftree_->Branch("z_lay_fit",&z_lay_fit_eff_);
  efftree_->Branch("ex_lay_fit",&ex_lay_fit_eff_);
  efftree_->Branch("ey_lay_fit",&ey_lay_fit_eff_);
  efftree_->Branch("ez_lay_fit",&ez_lay_fit_eff_);
  
  // end of position pull info

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

void TTreeValidation::collectSimTkTSVecMapInfo(const unsigned int mcTrackID, const TSVec& initTSs){
  simTkTSVecMap_[mcTrackID] = initTSs;
}

void TTreeValidation::collectSeedTkCFMapInfo(const unsigned int seedID, const TrackState& cfitStateHit0){
  seedTkCFMap_[seedID] = cfitStateHit0;
}

void TTreeValidation::collectSeedTkTSLayerPairVecMapInfo(const unsigned int seedID, const TSLayerPairVec& updatedStates){
  seedTkTSLayerPairVecMap_[seedID] = updatedStates;
}

void TTreeValidation::collectBranchingInfo(const unsigned int seedID, const unsigned int ilayer, const float nSigmaDeta, const float etaBinMinus, const unsigned int etaBinPlus, const float nSigmaDphi, const unsigned int phiBinMinus, const unsigned int phiBinPlus, const std::vector<unsigned int> & cand_hit_indices, const std::vector<unsigned int> branch_hit_indices){
  
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

void TTreeValidation::collectFitTkCFMapInfo(const unsigned int seedID, const TrackState& cfitStateHit0){
  fitTkCFMap_[seedID] = cfitStateHit0;
}

void TTreeValidation::collectFitTkTSLayerPairVecMapInfo(const unsigned int seedID, const TSLayerPairVec& updatedStates){
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

void TTreeValidation::fillSegmentTree(const BinInfoMap& segmentMap, const unsigned int evtID){
  std::lock_guard<std::mutex> locker(glock_);

  evtID_seg_ = evtID;
  for (unsigned int i = 0; i < Config::nLayers; i++) {
    layer_seg_ = i;
    for (unsigned int j = 0; j < Config::nEtaPart; j++) {
      etabin_seg_ = j;
      for (unsigned int k = 0; k < Config::nPhiPart; k++) {
	phibin_seg_ = k;
	nHits_seg_  = segmentMap[i][j][k].second;

	segtree_->Fill();
      }
    }
  }
}

void TTreeValidation::fillBranchTree(const unsigned int evtID)
{
  std::lock_guard<std::mutex> locker(glock_);
  
  evtID_br_  = evtID;
  for (TkToBVVMMIter seediter = seedToBranchValVecLayMapMap_.begin(); seediter != seedToBranchValVecLayMapMap_.end(); ++seediter){
    seedID_br_ = (*seediter).first;
    /*
    std::cout << std::endl;
    std::cout << "Seed ID: " << seedID_br_ << std::endl;
    */
    for (BVVLMiter layiter = (*seediter).second.begin(); layiter != (*seediter).second.end(); ++layiter){
      const auto& BranchValVec((*layiter).second);
      const unsigned int cands = BranchValVec.size();
 
      // clear vectors before filling
      candEtaPhiBins_.clear();
      candHits_.clear();
      candBranches_.clear();
      candnSigmaDeta_.clear();
      candnSigmaDphi_.clear();
        
      // totals
      std::vector<unsigned int> candEtaPhiBins(cands);
      std::vector<unsigned int> candHits(cands);
      std::vector<unsigned int> candBranches(cands);

      // unique hits, etaphibins, branches...
      std::unordered_map<unsigned int, bool> uniqueEtaPhiBins; // once a bin, branch, hit is used, set to true to count it only once. take size of map as "uniques"
      std::unordered_map<unsigned int, bool> uniqueHits;
      std::unordered_map<unsigned int, bool> uniqueBranches;
  
      // nSigmaDeta/phi vec
      std::vector<float> candnSigmaDeta(cands);
      std::vector<float> candnSigmaDphi(cands);

      for (unsigned int cand = 0; cand < cands; cand++){ // loop over input candidates at this layer for this seed
	const auto& BranchVal(BranchValVec[cand]); // grab the branch validation object
	
	/*	if ((*layiter).first == 9){
	  std::cout << "Cand: " << cand << std::endl;
	  std::cout << "pM: " << BranchVal.phiBinMinus << " pP: " << BranchVal.phiBinPlus << std::endl;
	  std::cout << "eM: " << BranchVal.etaBinMinus << " eP: " << BranchVal.etaBinPlus << std::endl;
	  }*/

	////////////////////////////////////
	//  EtaPhiBins Explored Counting  //
	////////////////////////////////////
	
	// set values for total etaphibins explored per candidate, also unique ones for all candidates
	if (BranchVal.phiBinPlus >= BranchVal.phiBinMinus){ // count the number of eta/phi bins explored if no phi wrapping
	  candEtaPhiBins[cand] = (BranchVal.etaBinPlus-BranchVal.etaBinMinus+1)*(BranchVal.phiBinPlus-BranchVal.phiBinMinus+1); // total etaphibins for this candidate

	  // no phi wrap to count uniques
	  for (unsigned int ibin = BranchVal.phiBinMinus; ibin <= BranchVal.phiBinPlus; ibin++){
	    uniqueEtaPhiBins[ibin] = true;
	  }
	}
	else{ // if phi wrapping, count need to do this obnoxious counting
	  candEtaPhiBins[cand] = (BranchVal.etaBinPlus-BranchVal.etaBinMinus+1)*(Config::nPhiPart-BranchVal.phiBinMinus+BranchVal.phiBinPlus+1); // total etaphibins for this candidate with phi wrapping

	  // use phi wrapping to count uniques
	  for (unsigned int ibin = BranchVal.phiBinMinus; ibin < Config::nPhiPart; ibin++){
	    uniqueEtaPhiBins[ibin] = true;
	  }
	  for (unsigned int ibin = 0; ibin <= BranchVal.phiBinPlus; ibin++){
	    uniqueEtaPhiBins[ibin] = true;
	  }
	}

	//////////////////////////////
	//  Hits Explored Counting  //
	//////////////////////////////

	candHits[cand] = BranchVal.cand_hit_indices.size(); // set value of nHits explored per input cand for this seed+layer
	for (auto&& cand_hit_idx : BranchVal.cand_hit_indices){ // save unique hits 
	  /*  if ((*layiter).first == 9){
	    std::cout << "layer9 hit: " << cand_hit_idx << std::endl;
	    }*/
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

      layer_  = (*layiter).first; // first index here is layer
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

void TTreeValidation::makeSimTkToRecoTksMaps(TrackVec& evt_seed_tracks, TrackVec& evt_build_tracks, TrackVec& evt_fit_tracks){
  std::lock_guard<std::mutex> locker(glock_);

  // set mcTkIDs... and sort by each (simTracks set in order by default!)
  mapSimTkToRecoTks(evt_seed_tracks,simToSeedMap_);
  mapSimTkToRecoTks(evt_build_tracks,simToBuildMap_);
  mapSimTkToRecoTks(evt_fit_tracks,simToFitMap_);
}

void TTreeValidation::mapSimTkToRecoTks(TrackVec& evt_tracks, TkToTkRefVecMap& simTkMap){
  for (auto&& track : evt_tracks){
    track.setMCTrackIDInfo();
    if (track.mcTrackID() != 999999){ // skip fakes, don't store them at all in sim map
      simTkMap[track.mcTrackID()].push_back(&track);
    }
  }

  for (auto&& simTkMatches : simTkMap){
    if (simTkMatches.second.size() < 2) { // no duplicates
      simTkMatches.second[0]->setMCDuplicateInfo(0,bool(false));
    }
    else{ // sort duplicates (ghosts) to keep best one --> most hits, lowest chi2
      std::sort(simTkMatches.second.begin(), simTkMatches.second.end(), sortByHitsChi2); 
      unsigned int duplicateID = 0;
      for (auto&& track : simTkMatches.second){
	track->setMCDuplicateInfo(duplicateID,bool(true));
	duplicateID++; // used in fake rate trees!
      } 
    }
  }
}

void TTreeValidation::makeSeedTkToRecoTkMaps(const TrackVec& evt_build_tracks, const TrackVec& evt_fit_tracks){
  std::lock_guard<std::mutex> locker(glock_); 

  // map seed to reco tracks --> seed track collection assumed to map to itself, unless we make some cuts
  mapSeedTkToRecoTk(evt_build_tracks,seedToBuildMap_);
  mapSeedTkToRecoTk(evt_fit_tracks,seedToFitMap_);
}

void TTreeValidation::mapSeedTkToRecoTk(const TrackVec& evt_tracks, TkToTkRefMap& seedTkMap){
  for (auto&& track : evt_tracks){
    seedTkMap[track.seedID()] = &track;
  }
}

void TTreeValidation::fillEffTree(const TrackVec& evt_sim_tracks, const unsigned int ievt){
  std::lock_guard<std::mutex> locker(glock_);

  for (auto&& simtrack : evt_sim_tracks){
    evtID_eff_ = ievt;
    mcID_eff_  = simtrack.mcTrackID();

    // generated values
    pt_mc_gen_eff_  = simtrack.pt(); 
    phi_mc_gen_eff_ = simtrack.momPhi();
    eta_mc_gen_eff_ = simtrack.momEta();
    nHits_mc_eff_   = simtrack.nHits();

    // for detector sim plots
    x_mc_reco_hit_eff_.clear();
    y_mc_reco_hit_eff_.clear();
    z_mc_reco_hit_eff_.clear();
    for (auto&& simhit : simtrack.hitsVector()){ // assume one hit per layer
      x_mc_reco_hit_eff_.push_back(simhit.x());
      y_mc_reco_hit_eff_.push_back(simhit.y());
      z_mc_reco_hit_eff_.push_back(simhit.z());
    }
    
    x_mc_gen_vrx_eff_ = simtrack.x();
    y_mc_gen_vrx_eff_ = simtrack.y();
    z_mc_gen_vrx_eff_ = simtrack.z();

    // clear vectors for position pulls
    layers_seed_eff_.clear();
    x_lay_seed_eff_.clear();
    y_lay_seed_eff_.clear();
    z_lay_seed_eff_.clear();
    ex_lay_seed_eff_.clear();
    ey_lay_seed_eff_.clear();
    ez_lay_seed_eff_.clear();

    layers_fit_eff_.clear();
    x_lay_fit_eff_.clear();
    y_lay_fit_eff_.clear();
    z_lay_fit_eff_.clear();
    ex_lay_fit_eff_.clear();
    ey_lay_fit_eff_.clear();
    ez_lay_fit_eff_.clear();

    // matched seed track
    if (simToSeedMap_.count(mcID_eff_)){ // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToSeedMap_[matched SimID][first element in vector]
      mcmask_seed_eff_ = 1; // quick logic for matched

      seedID_seed_eff_ = simToSeedMap_[mcID_eff_][0]->seedID(); 
      // use this to access correct sim track layer params
      const float layer = simToSeedMap_[mcID_eff_][0]->hitsVector().back().layer(); // last layer seed ended up on
      const SVector6 & initLayPrms = simTkTSVecMap_[mcID_eff_][layer].parameters; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps

      pt_mc_seed_eff_  = std::sqrt(getRad2(initLayPrms.At(3),initLayPrms.At(4)));
      phi_mc_seed_eff_ = getPhi(initLayPrms.At(3),initLayPrms.At(4));
      eta_mc_seed_eff_ = getEta(pt_mc_seed_eff_,initLayPrms.At(5));

      pt_seed_eff_   = simToSeedMap_[mcID_eff_][0]->pt(); 
      ept_seed_eff_  = simToSeedMap_[mcID_eff_][0]->ept();
      phi_seed_eff_  = simToSeedMap_[mcID_eff_][0]->momPhi(); 
      ephi_seed_eff_ = simToSeedMap_[mcID_eff_][0]->emomPhi();
      eta_seed_eff_  = simToSeedMap_[mcID_eff_][0]->momEta(); 
      eeta_seed_eff_ = simToSeedMap_[mcID_eff_][0]->emomEta();

      // Conformal fit stuff to match how it was done before
      const float cflayer = simToSeedMap_[mcID_eff_][0]->hitsVector().front().layer(); //  layer for which cf parameters are calculated with respect to
      const SVector6 & cfInitLayPrms = simTkTSVecMap_[mcID_eff_][cflayer].parameters; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps
      x_mc_cf_seed_eff_ = cfInitLayPrms.At(0);
      y_mc_cf_seed_eff_ = cfInitLayPrms.At(1);
      z_mc_cf_seed_eff_ = cfInitLayPrms.At(2);

      px_mc_cf_seed_eff_ = cfInitLayPrms.At(3);
      py_mc_cf_seed_eff_ = cfInitLayPrms.At(4);
      pz_mc_cf_seed_eff_ = cfInitLayPrms.At(5);

      pt_mc_cf_seed_eff_    = std::sqrt(getRad2(px_mc_cf_seed_eff_,py_mc_cf_seed_eff_));
      invpt_mc_cf_seed_eff_ = std::sqrt(getInvRad2(px_mc_cf_seed_eff_,py_mc_cf_seed_eff_));
      phi_mc_cf_seed_eff_   = getPhi(px_mc_cf_seed_eff_,py_mc_cf_seed_eff_);
      theta_mc_cf_seed_eff_ = getTheta(pt_mc_cf_seed_eff_,pz_mc_cf_seed_eff_);

      const TrackState & cfSeedTS = seedTkCFMap_[seedID_seed_eff_];

      x_cf_seed_eff_ = cfSeedTS.parameters.At(0);
      y_cf_seed_eff_ = cfSeedTS.parameters.At(1);
      z_cf_seed_eff_ = cfSeedTS.parameters.At(2);
      
      ex_cf_seed_eff_ = std::sqrt(cfSeedTS.errors.At(0,0));
      ey_cf_seed_eff_ = std::sqrt(cfSeedTS.errors.At(1,1));
      ez_cf_seed_eff_ = std::sqrt(cfSeedTS.errors.At(2,2));

      px_cf_seed_eff_ = cfSeedTS.parameters.At(3);
      py_cf_seed_eff_ = cfSeedTS.parameters.At(4);
      pz_cf_seed_eff_ = cfSeedTS.parameters.At(5);

      epx_cf_seed_eff_ = std::sqrt(cfSeedTS.errors.At(3,3));
      epy_cf_seed_eff_ = std::sqrt(cfSeedTS.errors.At(4,4));
      epz_cf_seed_eff_ = std::sqrt(cfSeedTS.errors.At(5,5));

      pt_cf_seed_eff_    = std::sqrt(getRad2(px_cf_seed_eff_,py_cf_seed_eff_));
      invpt_cf_seed_eff_ = std::sqrt(getInvRad2(px_cf_seed_eff_,py_cf_seed_eff_));
      phi_cf_seed_eff_   = getPhi(px_cf_seed_eff_,py_cf_seed_eff_);
      theta_cf_seed_eff_ = getTheta(pt_cf_seed_eff_,pz_cf_seed_eff_);

      ept_cf_seed_eff_    = std::sqrt(getRadErr2(px_cf_seed_eff_,py_cf_seed_eff_,epx_cf_seed_eff_,epy_cf_seed_eff_,cfSeedTS.errors(3,4)));
      einvpt_cf_seed_eff_ = std::sqrt(getInvRadErr2(px_cf_seed_eff_,py_cf_seed_eff_,epx_cf_seed_eff_,epy_cf_seed_eff_,cfSeedTS.errors(3,4)));
      ephi_cf_seed_eff_   = std::sqrt(getPhiErr2(px_cf_seed_eff_,py_cf_seed_eff_,epx_cf_seed_eff_,epy_cf_seed_eff_,cfSeedTS.errors(3,4)));
      etheta_cf_seed_eff_ = std::sqrt(getThetaErr2(cfSeedTS.parameters.At(3),cfSeedTS.parameters.At(4),cfSeedTS.parameters.At(5),cfSeedTS.errors.At(3,3),cfSeedTS.errors.At(4,4),cfSeedTS.errors(5,5),cfSeedTS.errors(3,4),cfSeedTS.errors(3,5),cfSeedTS.errors(4,5)));

      // position pull info
      const TSLayerPairVec & seedTSLayerPairVec = seedTkTSLayerPairVecMap_[seedID_seed_eff_];
      for (unsigned int ilay = 0; ilay < seedTSLayerPairVec.size(); ilay++){ // loop over layers present in pair vector
	layers_seed_eff_.push_back(seedTSLayerPairVec[ilay].first); // want to push back the ACTUAL layers stored, not the index of the loop!
	
	x_lay_seed_eff_.push_back(seedTSLayerPairVec[ilay].second.parameters.At(0)); // therefore, trackstate filled in sync with the layer it was saved on for the vector
	y_lay_seed_eff_.push_back(seedTSLayerPairVec[ilay].second.parameters.At(1));
	z_lay_seed_eff_.push_back(seedTSLayerPairVec[ilay].second.parameters.At(2));
	
	ex_lay_seed_eff_.push_back(std::sqrt(seedTSLayerPairVec[ilay].second.errors.At(0,0)));
	ey_lay_seed_eff_.push_back(std::sqrt(seedTSLayerPairVec[ilay].second.errors.At(1,1)));
	ez_lay_seed_eff_.push_back(std::sqrt(seedTSLayerPairVec[ilay].second.errors.At(2,2)));
      }

      // rest of mc info
      nHits_seed_eff_           = simToSeedMap_[mcID_eff_][0]->nHits();
      nHitsMatched_seed_eff_    = simToSeedMap_[mcID_eff_][0]->nHitsMatched();
      fracHitsMatched_seed_eff_ = float(nHitsMatched_seed_eff_) / float(nHits_seed_eff_);

      hitchi2_seed_eff_   = -10; //simToSeedMap_[mcID_eff_][0]->chi2(); // currently not being used
      helixchi2_seed_eff_ = computeHelixChi2(initLayPrms,simToSeedMap_[mcID_eff_][0]->parameters(),simToSeedMap_[mcID_eff_][0]->errors());

      duplmask_seed_eff_   = simToSeedMap_[mcID_eff_][0]->isDuplicate(); 
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

      x_mc_cf_seed_eff_ = -2000;
      y_mc_cf_seed_eff_ = -2000;
      z_mc_cf_seed_eff_ = -2000;

      px_mc_cf_seed_eff_ = -99;
      py_mc_cf_seed_eff_ = -99;
      pz_mc_cf_seed_eff_ = -99;

      pt_mc_cf_seed_eff_    = -99;
      invpt_mc_cf_seed_eff_ = -99;
      phi_mc_cf_seed_eff_   = -99;
      theta_mc_cf_seed_eff_ = -99;

      x_cf_seed_eff_ = -2000;
      y_cf_seed_eff_ = -2000;
      z_cf_seed_eff_ = -2000;
      
      ex_cf_seed_eff_ = -2000;
      ey_cf_seed_eff_ = -2000;
      ez_cf_seed_eff_ = -2000;

      px_cf_seed_eff_ = -99;
      py_cf_seed_eff_ = -99;
      pz_cf_seed_eff_ = -99;

      epx_cf_seed_eff_ = -99;
      epy_cf_seed_eff_ = -99;
      epz_cf_seed_eff_ = -99;

      pt_cf_seed_eff_    = -99;
      invpt_cf_seed_eff_ = -99;
      phi_cf_seed_eff_   = -99;
      theta_cf_seed_eff_ = -99;

      ept_cf_seed_eff_    = -99;
      einvpt_cf_seed_eff_ = -99;
      ephi_cf_seed_eff_   = -99;
      etheta_cf_seed_eff_ = -99;

      layers_seed_eff_.push_back(-1); // mask per layer
      
      x_lay_seed_eff_.push_back(-2000);
      y_lay_seed_eff_.push_back(-2000);
      z_lay_seed_eff_.push_back(-2000);
      
      ex_lay_seed_eff_.push_back(-2000);
      ey_lay_seed_eff_.push_back(-2000);
      ez_lay_seed_eff_.push_back(-2000);

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
      mcmask_build_eff_ = 1; // quick logic for matched

      seedID_build_eff_ = simToBuildMap_[mcID_eff_][0]->seedID();
      // use this to access correct sim track layer params
      const float layer = simToBuildMap_[mcID_eff_][0]->hitsVector().back().layer(); // last layer build track ended up on
      const SVector6 & initLayPrms = simTkTSVecMap_[mcID_eff_][layer].parameters; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps

      pt_mc_build_eff_  = std::sqrt(getRad2(initLayPrms.At(3),initLayPrms.At(4)));
      phi_mc_build_eff_ = getPhi(initLayPrms.At(3),initLayPrms.At(4));
      eta_mc_build_eff_ = getEta(pt_mc_build_eff_,initLayPrms.At(5));

      pt_build_eff_   = simToBuildMap_[mcID_eff_][0]->pt(); 
      ept_build_eff_  = simToBuildMap_[mcID_eff_][0]->ept();
      phi_build_eff_  = simToBuildMap_[mcID_eff_][0]->momPhi(); 
      ephi_build_eff_ = simToBuildMap_[mcID_eff_][0]->emomPhi();
      eta_build_eff_  = simToBuildMap_[mcID_eff_][0]->momEta(); 
      eeta_build_eff_ = simToBuildMap_[mcID_eff_][0]->emomEta();
   
      nHits_build_eff_           = simToBuildMap_[mcID_eff_][0]->nHits();
      nHitsMatched_build_eff_    = simToBuildMap_[mcID_eff_][0]->nHitsMatched();
      fracHitsMatched_build_eff_ = float(nHitsMatched_build_eff_) / float(nHits_build_eff_);

      hitchi2_build_eff_   = -10; //simToBuildMap_[mcID_eff_][0]->chi2(); // currently not being used
      helixchi2_build_eff_ = computeHelixChi2(initLayPrms,simToBuildMap_[mcID_eff_][0]->parameters(),simToBuildMap_[mcID_eff_][0]->errors());

      duplmask_build_eff_   = simToBuildMap_[mcID_eff_][0]->isDuplicate(); 
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
      mcmask_fit_eff_ = 1; // quick logic for matched

      seedID_fit_eff_ = simToFitMap_[mcID_eff_][0]->seedID();
      // use this to access correct sim track layer params
      const float layer = simToFitMap_[mcID_eff_][0]->hitsVector().back().layer(); // last layer fit track ended up on
      const SVector6 & initLayPrms = simTkTSVecMap_[mcID_eff_][layer].parameters; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps

      pt_mc_fit_eff_  = std::sqrt(getRad2(initLayPrms.At(3),initLayPrms.At(4)));
      phi_mc_fit_eff_ = getPhi(initLayPrms.At(3),initLayPrms.At(4));
      eta_mc_fit_eff_ = getEta(pt_mc_fit_eff_,initLayPrms.At(5));

      pt_fit_eff_   = simToFitMap_[mcID_eff_][0]->pt(); 
      ept_fit_eff_  = simToFitMap_[mcID_eff_][0]->ept();
      phi_fit_eff_  = simToFitMap_[mcID_eff_][0]->momPhi(); 
      ephi_fit_eff_ = simToFitMap_[mcID_eff_][0]->emomPhi();
      eta_fit_eff_  = simToFitMap_[mcID_eff_][0]->momEta(); 
      eeta_fit_eff_ = simToFitMap_[mcID_eff_][0]->emomEta();

      // Conformal utils same info for plots
      const float cflayer = simToFitMap_[mcID_eff_][0]->hitsVector().front().layer();
      const SVector6 & cfInitLayPrms = simTkTSVecMap_[mcID_eff_][cflayer].parameters; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps
      x_mc_cf_fit_eff_ = cfInitLayPrms.At(0);
      y_mc_cf_fit_eff_ = cfInitLayPrms.At(1);
      z_mc_cf_fit_eff_ = cfInitLayPrms.At(2);

      px_mc_cf_fit_eff_ = cfInitLayPrms.At(3);
      py_mc_cf_fit_eff_ = cfInitLayPrms.At(4);
      pz_mc_cf_fit_eff_ = cfInitLayPrms.At(5);

      pt_mc_cf_fit_eff_    = std::sqrt(getRad2(px_mc_cf_fit_eff_,py_mc_cf_fit_eff_));
      invpt_mc_cf_fit_eff_ = std::sqrt(getInvRad2(px_mc_cf_fit_eff_,py_mc_cf_fit_eff_));
      phi_mc_cf_fit_eff_   = getPhi(px_mc_cf_fit_eff_,py_mc_cf_fit_eff_);
      theta_mc_cf_fit_eff_ = getTheta(pt_mc_cf_fit_eff_,pz_mc_cf_fit_eff_);

      const TrackState & cfFitTS = fitTkCFMap_[seedID_fit_eff_];

      x_cf_fit_eff_ = cfFitTS.parameters.At(0);
      y_cf_fit_eff_ = cfFitTS.parameters.At(1);
      z_cf_fit_eff_ = cfFitTS.parameters.At(2);
      
      ex_cf_fit_eff_ = std::sqrt(cfFitTS.errors.At(0,0));
      ey_cf_fit_eff_ = std::sqrt(cfFitTS.errors.At(1,1));
      ez_cf_fit_eff_ = std::sqrt(cfFitTS.errors.At(2,2));

      px_cf_fit_eff_ = cfFitTS.parameters.At(3);
      py_cf_fit_eff_ = cfFitTS.parameters.At(4);
      pz_cf_fit_eff_ = cfFitTS.parameters.At(5);

      epx_cf_fit_eff_ = std::sqrt(cfFitTS.errors.At(3,3));
      epy_cf_fit_eff_ = std::sqrt(cfFitTS.errors.At(4,4));
      epz_cf_fit_eff_ = std::sqrt(cfFitTS.errors.At(5,5));

      pt_cf_fit_eff_    = std::sqrt(getRad2(px_cf_fit_eff_,py_cf_fit_eff_));
      invpt_cf_fit_eff_ = std::sqrt(getInvRad2(px_cf_fit_eff_,py_cf_fit_eff_));
      phi_cf_fit_eff_   = getPhi(px_cf_fit_eff_,py_cf_fit_eff_);
      theta_cf_fit_eff_ = getTheta(pt_cf_fit_eff_,pz_cf_fit_eff_);

      ept_cf_fit_eff_    = std::sqrt(getRadErr2(cfFitTS.parameters.At(3),cfFitTS.parameters.At(4),cfFitTS.errors.At(3,3),cfFitTS.errors.At(4,4),cfFitTS.errors(3,4)));
      einvpt_cf_fit_eff_ = std::sqrt(getInvRadErr2(cfFitTS.parameters.At(3),cfFitTS.parameters.At(4),cfFitTS.errors.At(3,3),cfFitTS.errors.At(4,4),cfFitTS.errors(3,4)));
      ephi_cf_fit_eff_   = std::sqrt(getPhiErr2(cfFitTS.parameters.At(3),cfFitTS.parameters.At(4),cfFitTS.errors.At(3,3),cfFitTS.errors.At(4,4),cfFitTS.errors(3,4)));
      etheta_cf_fit_eff_ = std::sqrt(getThetaErr2(cfFitTS.parameters.At(3),cfFitTS.parameters.At(4),cfFitTS.parameters.At(5),cfFitTS.errors.At(3,3),cfFitTS.errors.At(4,4),cfFitTS.errors(5,5),cfFitTS.errors(3,4),cfFitTS.errors(3,5),cfFitTS.errors(4,5)));
      
      // position pull info
      const TSLayerPairVec & fitTSLayerPairVec = fitTkTSLayerPairVecMap_[seedID_fit_eff_];
      for (unsigned int ilay = 0; ilay < fitTSLayerPairVec.size(); ilay++){
	layers_fit_eff_.push_back(fitTSLayerPairVec[ilay].first); // want ACTUAL layers fit landed on, not the index of pair vec
	
	x_lay_fit_eff_.push_back(fitTSLayerPairVec[ilay].second.parameters.At(0)); // TS in sync with layer it landed on
	y_lay_fit_eff_.push_back(fitTSLayerPairVec[ilay].second.parameters.At(1));
	z_lay_fit_eff_.push_back(fitTSLayerPairVec[ilay].second.parameters.At(2));
	
	ex_lay_fit_eff_.push_back(std::sqrt(fitTSLayerPairVec[ilay].second.errors.At(0,0)));
	ey_lay_fit_eff_.push_back(std::sqrt(fitTSLayerPairVec[ilay].second.errors.At(1,1)));
	ez_lay_fit_eff_.push_back(std::sqrt(fitTSLayerPairVec[ilay].second.errors.At(2,2)));
      }

      // rest of mc info
      nHits_fit_eff_           = simToFitMap_[mcID_eff_][0]->nHits();
      nHitsMatched_fit_eff_    = simToFitMap_[mcID_eff_][0]->nHitsMatched();
      fracHitsMatched_fit_eff_ = float(nHitsMatched_fit_eff_) / float(nHits_fit_eff_);

      hitchi2_fit_eff_   = -10; //simToFitMap_[mcID_eff_][0]->chi2(); // currently not being used
      helixchi2_fit_eff_ = computeHelixChi2(initLayPrms,simToFitMap_[mcID_eff_][0]->parameters(),simToFitMap_[mcID_eff_][0]->errors());

      duplmask_fit_eff_   = simToFitMap_[mcID_eff_][0]->isDuplicate(); 
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

      x_mc_cf_fit_eff_ = -2000;
      y_mc_cf_fit_eff_ = -2000;
      z_mc_cf_fit_eff_ = -2000;

      px_mc_cf_fit_eff_ = -99;
      py_mc_cf_fit_eff_ = -99;
      pz_mc_cf_fit_eff_ = -99;

      pt_mc_cf_fit_eff_    = -99;
      invpt_mc_cf_fit_eff_ = -99;
      phi_mc_cf_fit_eff_   = -99;
      theta_mc_cf_fit_eff_ = -99;

      x_cf_fit_eff_ = -2000;
      y_cf_fit_eff_ = -2000;
      z_cf_fit_eff_ = -2000;
      
      ex_cf_fit_eff_ = -2000;
      ey_cf_fit_eff_ = -2000;
      ez_cf_fit_eff_ = -2000;

      px_cf_fit_eff_ = -99;
      py_cf_fit_eff_ = -99;
      pz_cf_fit_eff_ = -99;

      epx_cf_fit_eff_ = -99;
      epy_cf_fit_eff_ = -99;
      epz_cf_fit_eff_ = -99;

      pt_cf_fit_eff_    = -99;
      invpt_cf_fit_eff_ = -99;
      phi_cf_fit_eff_   = -99;
      theta_cf_fit_eff_ = -99;

      ept_cf_fit_eff_    = -99;
      einvpt_cf_fit_eff_ = -99;
      ephi_cf_fit_eff_   = -99;
      etheta_cf_fit_eff_ = -99;

      layers_fit_eff_.push_back(-1); // mask per layer
      
      x_lay_fit_eff_.push_back(-2000);
      y_lay_fit_eff_.push_back(-2000);
      z_lay_fit_eff_.push_back(-2000);
      
      ex_lay_fit_eff_.push_back(-2000);
      ey_lay_fit_eff_.push_back(-2000);
      ez_lay_fit_eff_.push_back(-2000);

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

void TTreeValidation::fillFakeRateTree(const TrackVec& evt_seed_tracks, const unsigned int ievt){
  std::lock_guard<std::mutex> locker(glock_);
  
  for (auto&& seedtrack : evt_seed_tracks){
    evtID_FR_  = ievt;
    seedID_FR_ = seedtrack.seedID();

    // seed info
    seedmask_seed_FR_ = 1; // automatically set to 1, because at the moment no cuts on seeds after conformal+KF fit.  seed triplets filetered by RZ chi2 before fitting. 

    pt_seed_FR_   = seedtrack.pt();
    ept_seed_FR_  = seedtrack.ept();
    phi_seed_FR_  = seedtrack.momPhi();
    ephi_seed_FR_ = seedtrack.emomPhi();
    eta_seed_FR_  = seedtrack.momEta();
    eeta_seed_FR_ = seedtrack.emomEta();

    nHits_seed_FR_           = seedtrack.nHits();
    nHitsMatched_seed_FR_    = seedtrack.nHitsMatched();
    fracHitsMatched_seed_FR_ = float(nHitsMatched_seed_FR_) / float(nHits_seed_FR_);

    hitchi2_seed_FR_ = -10; // seedtrack.chi2(); --> not currently used

    // sim info for seed track
    mcID_seed_FR_ = seedtrack.mcTrackID();
    if (mcID_seed_FR_ != 999999){ // store sim info at that final layer!!! --> gen info stored only in eff tree
      mcmask_seed_FR_ = 1; // matched track to sim

      const float layer = seedtrack.hitsVector().back().layer();
      const SVector6 & initLayPrms = simTkTSVecMap_[mcID_seed_FR_][layer].parameters; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps
      
      pt_mc_seed_FR_    = std::sqrt(getRad2(initLayPrms.At(3),initLayPrms.At(4))); // again, pt of mc truth at the layer the seed ends
      phi_mc_seed_FR_   = getPhi(initLayPrms.At(3),initLayPrms.At(4));
      eta_mc_seed_FR_   = getEta(pt_mc_seed_FR_,initLayPrms.At(5));
      nHits_mc_seed_FR_ = simTkTSVecMap_[mcID_seed_FR_].size();

      helixchi2_seed_FR_ = computeHelixChi2(initLayPrms,seedtrack.parameters(),seedtrack.errors());

      duplmask_seed_FR_   = seedtrack.isDuplicate();
      iTkMatches_seed_FR_ = seedtrack.duplicateID(); // ith duplicate seed track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"      
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

      pt_build_FR_   = seedToBuildMap_[seedID_FR_]->pt();
      ept_build_FR_  = seedToBuildMap_[seedID_FR_]->ept();
      phi_build_FR_  = seedToBuildMap_[seedID_FR_]->momPhi();
      ephi_build_FR_ = seedToBuildMap_[seedID_FR_]->emomPhi();
      eta_build_FR_  = seedToBuildMap_[seedID_FR_]->momEta();
      eeta_build_FR_ = seedToBuildMap_[seedID_FR_]->emomEta();

      nHits_build_FR_           = seedToBuildMap_[seedID_FR_]->nHits();
      nHitsMatched_build_FR_    = seedToBuildMap_[seedID_FR_]->nHitsMatched();
      fracHitsMatched_build_FR_ = float(nHitsMatched_build_FR_) / float(nHits_build_FR_);

      hitchi2_build_FR_ = seedToBuildMap_[seedID_FR_]->chi2();

      // sim info for build track
      mcID_build_FR_  = seedToBuildMap_[seedID_FR_]->mcTrackID();
      if (mcID_build_FR_ != 999999){ // build track matched to seed and sim 
	mcmask_build_FR_ = 1; // matched track to sim

	const float layer = seedToBuildMap_[seedID_FR_]->hitsVector().back().layer(); // last layer building ended up on
	const SVector6 & initLayPrms = simTkTSVecMap_[mcID_build_FR_][layer].parameters; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps
      
	pt_mc_build_FR_    = std::sqrt(getRad2(initLayPrms.At(3),initLayPrms.At(4))); // again, pt of mc truth at the layer the seed ends
	phi_mc_build_FR_   = getPhi(initLayPrms.At(3),initLayPrms.At(4));
	eta_mc_build_FR_   = getEta(pt_mc_build_FR_,initLayPrms.At(5));
	nHits_mc_build_FR_ = simTkTSVecMap_[mcID_build_FR_].size(); // size of this vector is same size as nHits 

	helixchi2_build_FR_ = computeHelixChi2(initLayPrms,seedToBuildMap_[seedID_FR_]->parameters(),seedToBuildMap_[seedID_FR_]->errors());

	duplmask_build_FR_   = seedToBuildMap_[seedID_FR_]->isDuplicate();
	iTkMatches_build_FR_ = seedToBuildMap_[seedID_FR_]->duplicateID(); // ith duplicate build track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"
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

      pt_fit_FR_   = seedToFitMap_[seedID_FR_]->pt();
      ept_fit_FR_  = seedToFitMap_[seedID_FR_]->ept();
      phi_fit_FR_  = seedToFitMap_[seedID_FR_]->momPhi();
      ephi_fit_FR_ = seedToFitMap_[seedID_FR_]->emomPhi();
      eta_fit_FR_  = seedToFitMap_[seedID_FR_]->momEta();
      eeta_fit_FR_ = seedToFitMap_[seedID_FR_]->emomEta();

      nHits_fit_FR_           = seedToFitMap_[seedID_FR_]->nHits();
      nHitsMatched_fit_FR_    = seedToFitMap_[seedID_FR_]->nHitsMatched();
      fracHitsMatched_fit_FR_ = float(nHitsMatched_fit_FR_) / float(nHits_fit_FR_);

      hitchi2_fit_FR_ = seedToFitMap_[seedID_FR_]->chi2();

      // sim info for fit track
      mcID_fit_FR_  = seedToFitMap_[seedID_FR_]->mcTrackID();
      if (mcID_fit_FR_ != 999999){ // fit track matched to seed and sim 
	mcmask_fit_FR_ = 1; // matched track to sim

	const float layer = seedToFitMap_[seedID_FR_]->hitsVector().back().layer(); // last layer fiting ended up on
	const SVector6 & initLayPrms = simTkTSVecMap_[mcID_fit_FR_][layer].parameters; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps
   
	pt_mc_fit_FR_    = std::sqrt(getRad2(initLayPrms.At(3),initLayPrms.At(4))); // again, pt of mc truth at the layer the seed ends
	phi_mc_fit_FR_   = getPhi(initLayPrms.At(3),initLayPrms.At(4));
	eta_mc_fit_FR_   = getEta(pt_mc_fit_FR_,initLayPrms.At(5));
	nHits_mc_fit_FR_ = simTkTSVecMap_[mcID_fit_FR_].size();

	helixchi2_fit_FR_ = computeHelixChi2(initLayPrms,seedToFitMap_[seedID_FR_]->parameters(),seedToFitMap_[seedID_FR_]->errors());

	duplmask_fit_FR_   = seedToFitMap_[seedID_FR_]->isDuplicate();
	iTkMatches_fit_FR_ = seedToFitMap_[seedID_FR_]->duplicateID(); // ith duplicate fit track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"
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
  maxCand_ = Config::maxCand;
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
