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

  if (!Config::super_debug) { // regular validation
    initializeSeedInfoTree();
    initializeSeedTree();
    initializeSegmentTree();
    initializeBranchTree();
    initializeEfficiencyTree();
    initializeFakeRateTree();  
    initializeGeometryTree();
    initializeConformalTree();
  }
  else {
    initializeDebugTree();
  }
  initializeConfigTree();
}

void TTreeValidation::initializeDebugTree(){
  nlayers_debug_ = Config::nLayers;

  // debug tree(lots and lots of variables)
  debugtree_ = new TTree("debugtree","debugtree");

  debugtree_->Branch("nlayers",&nlayers_debug_,"nlayers_debug_/I");

  debugtree_->Branch("recocharge",&recocharge_debug_,"recocharge/I");
  debugtree_->Branch("mccharge",&mccharge_debug_,"mccharge/I");

  debugtree_->Branch("event",&event_debug_,"event/I");
  debugtree_->Branch("nHits",&nHits_debug_,"nHits/I");

  debugtree_->Branch("pt_gen",&pt_gen_debug_,"pt_gen/F");
  debugtree_->Branch("phi_gen",&phi_gen_debug_,"phi_gen/F");
  debugtree_->Branch("eta_gen",&eta_gen_debug_,"eta_gen/F");

  debugtree_->Branch("x_gen",&x_gen_debug_,"x_gen/F");
  debugtree_->Branch("y_gen",&y_gen_debug_,"y_gen/F");
  debugtree_->Branch("z_gen",&z_gen_debug_,"z_gen/F");

  debugtree_->Branch("layer_mc",&layer_mc_debug_,"layer_mc[nlayers_debug_]/I");

  debugtree_->Branch("layer_chi2",&layer_chi2_debug_,"layer_chi2[nlayers_debug_]/I");
  debugtree_->Branch("chi2",&chi2_debug_,"chi2[nlayers_debug_]/F");

  // MC
  debugtree_->Branch("x_hit",&x_hit_debug_,"x_hit[nlayers_debug_]/F");
  debugtree_->Branch("y_hit",&y_hit_debug_,"y_hit[nlayers_debug_]/F");
  debugtree_->Branch("z_hit",&z_hit_debug_,"z_hit[nlayers_debug_]/F");
  debugtree_->Branch("exx_hit",&exx_hit_debug_,"exx_hit[nlayers_debug_]/F");
  debugtree_->Branch("eyy_hit",&eyy_hit_debug_,"eyy_hit[nlayers_debug_]/F");
  debugtree_->Branch("ezz_hit",&ezz_hit_debug_,"ezz_hit[nlayers_debug_]/F");

  debugtree_->Branch("x_mc",&x_mc_debug_,"x_mc[nlayers_debug_]/F");
  debugtree_->Branch("y_mc",&y_mc_debug_,"y_mc[nlayers_debug_]/F");
  debugtree_->Branch("z_mc",&z_mc_debug_,"z_mc[nlayers_debug_]/F");
  debugtree_->Branch("px_mc",&px_mc_debug_,"px_mc[nlayers_debug_]/F");
  debugtree_->Branch("py_mc",&py_mc_debug_,"py_mc[nlayers_debug_]/F");
  debugtree_->Branch("pz_mc",&pz_mc_debug_,"pz_mc[nlayers_debug_]/F");

  debugtree_->Branch("pt_mc",&pt_mc_debug_,"pt_mc[nlayers_debug_]/F");
  debugtree_->Branch("phi_mc",&phi_mc_debug_,"phi_mc[nlayers_debug_]/F");
  debugtree_->Branch("eta_mc",&eta_mc_debug_,"eta_mc[nlayers_debug_]/F");
  debugtree_->Branch("invpt_mc",&invpt_mc_debug_,"invpt_mc[nlayers_debug_]/F");
  debugtree_->Branch("theta_mc",&theta_mc_debug_,"theta_mc[nlayers_debug_]/F");

  // conformal
  debugtree_->Branch("x_cf",&x_cf_debug_,"x_cf/F");
  debugtree_->Branch("y_cf",&y_cf_debug_,"y_cf/F");
  debugtree_->Branch("z_cf",&z_cf_debug_,"z_cf/F");
  debugtree_->Branch("exx_cf",&exx_cf_debug_,"exx_cf/F");
  debugtree_->Branch("eyy_cf",&eyy_cf_debug_,"eyy_cf/F");
  debugtree_->Branch("ezz_cf",&ezz_cf_debug_,"ezz_cf/F");

  debugtree_->Branch("px_cf",&px_cf_debug_,"px_cf/F");
  debugtree_->Branch("py_cf",&py_cf_debug_,"py_cf/F");
  debugtree_->Branch("pz_cf",&pz_cf_debug_,"pz_cf/F");
  debugtree_->Branch("epxpx_cf",&epxpx_cf_debug_,"epxpx_cf/F");
  debugtree_->Branch("epypy_cf",&epypy_cf_debug_,"epypy_cf/F");
  debugtree_->Branch("epzpz_cf",&epzpz_cf_debug_,"epzpz_cf/F");

  debugtree_->Branch("pt_cf",&pt_cf_debug_,"pt_cf/F");
  debugtree_->Branch("phi_cf",&phi_cf_debug_,"phi_cf/F");
  debugtree_->Branch("eta_cf",&eta_cf_debug_,"eta_cf/F");
  debugtree_->Branch("ept_cf",&ept_cf_debug_,"ept_cf/F");
  debugtree_->Branch("ephi_cf",&ephi_cf_debug_,"ephi_cf/F");
  debugtree_->Branch("eeta_cf",&eeta_cf_debug_,"eeta_cf/F");

  debugtree_->Branch("invpt_cf",&invpt_cf_debug_,"invpt_cf/F");
  debugtree_->Branch("einvpt_cf",&einvpt_cf_debug_,"einvpt_cf/F");
  debugtree_->Branch("theta_cf",&theta_cf_debug_,"theta_cf/F");
  debugtree_->Branch("etheta_cf",&etheta_cf_debug_,"etheta_cf/F");

  // prop
  debugtree_->Branch("layer_prop",&layer_prop_debug_,"layer_prop[nlayers_debug_]/I");

  debugtree_->Branch("x_prop",&x_prop_debug_,"x_prop[nlayers_debug_]/F");
  debugtree_->Branch("y_prop",&y_prop_debug_,"y_prop[nlayers_debug_]/F");
  debugtree_->Branch("z_prop",&z_prop_debug_,"z_prop[nlayers_debug_]/F");
  debugtree_->Branch("exx_prop",&exx_prop_debug_,"exx_prop[nlayers_debug_]/F");
  debugtree_->Branch("eyy_prop",&eyy_prop_debug_,"eyy_prop[nlayers_debug_]/F");
  debugtree_->Branch("ezz_prop",&ezz_prop_debug_,"ezz_prop[nlayers_debug_]/F");

  debugtree_->Branch("px_prop",&px_prop_debug_,"px_prop[nlayers_debug_]/F");
  debugtree_->Branch("py_prop",&py_prop_debug_,"py_prop[nlayers_debug_]/F");
  debugtree_->Branch("pz_prop",&pz_prop_debug_,"pz_prop[nlayers_debug_]/F");
  debugtree_->Branch("epxpx_prop",&epxpx_prop_debug_,"epxpx_prop[nlayers_debug_]/F");
  debugtree_->Branch("epypy_prop",&epypy_prop_debug_,"epypy_prop[nlayers_debug_]/F");
  debugtree_->Branch("epzpz_prop",&epzpz_prop_debug_,"epzpz_prop[nlayers_debug_]/F");

  debugtree_->Branch("pt_prop",&pt_prop_debug_,"pt_prop[nlayers_debug_]/F");
  debugtree_->Branch("phi_prop",&phi_prop_debug_,"phi_prop[nlayers_debug_]/F");
  debugtree_->Branch("eta_prop",&eta_prop_debug_,"eta_prop[nlayers_debug_]/F");
  debugtree_->Branch("ept_prop",&ept_prop_debug_,"ept_prop[nlayers_debug_]/F");
  debugtree_->Branch("ephi_prop",&ephi_prop_debug_,"ephi_prop[nlayers_debug_]/F");
  debugtree_->Branch("eeta_prop",&eeta_prop_debug_,"eeta_prop[nlayers_debug_]/F");

  debugtree_->Branch("invpt_prop",&invpt_prop_debug_,"invpt_prop[nlayers_debug_]/F");
  debugtree_->Branch("einvpt_prop",&einvpt_prop_debug_,"einvpt_prop[nlayers_debug_]/F");
  debugtree_->Branch("theta_prop",&theta_prop_debug_,"theta_prop[nlayers_debug_]/F");
  debugtree_->Branch("etheta_prop",&etheta_prop_debug_,"etheta_prop[nlayers_debug_]/F");

  //update
  debugtree_->Branch("layer_up",&layer_up_debug_,"layer_up[nlayers_debug_]/I");

  debugtree_->Branch("x_up",&x_up_debug_,"x_up[nlayers_debug_]/F");
  debugtree_->Branch("y_up",&y_up_debug_,"y_up[nlayers_debug_]/F");
  debugtree_->Branch("z_up",&z_up_debug_,"z_up[nlayers_debug_]/F");
  debugtree_->Branch("exx_up",&exx_up_debug_,"exx_up[nlayers_debug_]/F");
  debugtree_->Branch("eyy_up",&eyy_up_debug_,"eyy_up[nlayers_debug_]/F");
  debugtree_->Branch("ezz_up",&ezz_up_debug_,"ezz_up[nlayers_debug_]/F");

  debugtree_->Branch("px_up",&px_up_debug_,"px_up[nlayers_debug_]/F");
  debugtree_->Branch("py_up",&py_up_debug_,"py_up[nlayers_debug_]/F");
  debugtree_->Branch("pz_up",&pz_up_debug_,"pz_up[nlayers_debug_]/F");
  debugtree_->Branch("epxpx_up",&epxpx_up_debug_,"epxpx_up[nlayers_debug_]/F");
  debugtree_->Branch("epypy_up",&epypy_up_debug_,"epypy_up[nlayers_debug_]/F");
  debugtree_->Branch("epzpz_up",&epzpz_up_debug_,"epzpz_up[nlayers_debug_]/F");

  debugtree_->Branch("pt_up",&pt_up_debug_,"pt_up[nlayers_debug_]/F");
  debugtree_->Branch("phi_up",&phi_up_debug_,"phi_up[nlayers_debug_]/F");
  debugtree_->Branch("eta_up",&eta_up_debug_,"eta_up[nlayers_debug_]/F");
  debugtree_->Branch("ept_up",&ept_up_debug_,"ept_up[nlayers_debug_]/F");
  debugtree_->Branch("ephi_up",&ephi_up_debug_,"ephi_up[nlayers_debug_]/F");
  debugtree_->Branch("eeta_up",&eeta_up_debug_,"eeta_up[nlayers_debug_]/F");

  debugtree_->Branch("invpt_up",&invpt_up_debug_,"invpt_up[nlayers_debug_]/F");
  debugtree_->Branch("einvpt_up",&einvpt_up_debug_,"einvpt_up[nlayers_debug_]/F");
  debugtree_->Branch("theta_up",&theta_up_debug_,"theta_up[nlayers_debug_]/F");
  debugtree_->Branch("etheta_up",&etheta_up_debug_,"etheta_up[nlayers_debug_]/F");

  debugtree_->Branch("etabin_hit",&ebhit_debug_,"etabin_hit[nlayers_debug_]/I");
  debugtree_->Branch("etabinplus",&ebp_debug_,"etabinplus[nlayers_debug_]/I");
  debugtree_->Branch("etabinminus",&ebm_debug_,"etabinminus[nlayers_debug_]/I");
  debugtree_->Branch("phibin_hit",&pbhit_debug_,"phibin_hit[nlayers_debug_]/I");
  debugtree_->Branch("phibinplus",&pbp_debug_,"phibinplus[nlayers_debug_]/I");
  debugtree_->Branch("phibinminus",pbm_debug_,"phibinminus[nlayers_debug_]/I");
}

void TTreeValidation::initializeSeedInfoTree(){
  // seed validation
  seedinfotree_ = new TTree("seedinfotree","seedinfotree");
  seedinfotree_->Branch("evtID",&evtID_seedinfo_);
  seedinfotree_->Branch("mcID",&mcID_seedinfo_);
  seedinfotree_->Branch("pt_gen",&pt_gen_seedinfo_);
  seedinfotree_->Branch("eta_gen",&eta_gen_seedinfo_);
  seedinfotree_->Branch("phi_gen",&phi_gen_seedinfo_);
  seedinfotree_->Branch("a",&a);
  seedinfotree_->Branch("b",&b);
  seedinfotree_->Branch("r",&r);
  seedinfotree_->Branch("d0",&d0);
  seedinfotree_->Branch("pass",&pass);
}

void TTreeValidation::initializeSeedTree(){
  // seed validation
  seedtree_ = new TTree("seedtree","seedtree");
  seedtree_->Branch("evtID",&evtID_seed_);
  seedtree_->Branch("nTkAll",&nTkAll_);
  seedtree_->Branch("nTkAllMC",&nTkAllMC_);
  seedtree_->Branch("nTkCut",&nTkCut_);
  seedtree_->Branch("nTkCutMC",&nTkCutMC_);
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

void TTreeValidation::collectPropTSLayerVecInfo(int layer, const TrackState& propTS){
  propTSLayerPairVec_.push_back(std::make_pair(layer,propTS)); 
}

void TTreeValidation::collectChi2LayerVecInfo(int layer, float chi2){
  chi2LayerPairVec_.push_back(std::make_pair(layer,chi2)); 
}

void TTreeValidation::collectUpTSLayerVecInfo(int layer, const TrackState& upTS){
  upTSLayerPairVec_.push_back(std::make_pair(layer,upTS)); 
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

void TTreeValidation::resetDebugVectors(){
  std::lock_guard<std::mutex> locker(glock_);
  
  propTSLayerPairVec_.clear(); // used exclusively for debugtree (prop states)
  chi2LayerPairVec_.clear();   // used exclusively for debugtree 
  upTSLayerPairVec_.clear();   // used exclusively for debugtree (updated states)
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

void TTreeValidation::resetDebugTreeArrays(){
  for (int i = 0; i < Config::nLayers; i++){
    // reset MC info
    layer_mc_debug_[i]=-99;
    x_hit_debug_[i]=-99;     y_hit_debug_[i]=-99;     z_hit_debug_[i]=-99; 
    exx_hit_debug_[i]=-99;   eyy_hit_debug_[i]=-99;   ezz_hit_debug_[i]=-99;
    x_mc_debug_[i]=-99;     y_mc_debug_[i]=-99;     z_mc_debug_[i]=-99; 
    px_mc_debug_[i]=-99;    py_mc_debug_[i]=-99;    pz_mc_debug_[i]=-99;
    pt_mc_debug_[i]=-99;    phi_mc_debug_[i]=-99;   eta_mc_debug_[i]=-99;
    invpt_mc_debug_[i]=-99; theta_mc_debug_[i]=-99;

    // reset prop info
    layer_prop_debug_[i]=-99;
    x_prop_debug_[i]=-99;      y_prop_debug_[i]=-99;      z_prop_debug_[i]=-99;
    exx_prop_debug_[i]=-99;    eyy_prop_debug_[i]=-99;    ezz_prop_debug_[i]=-99;
    px_prop_debug_[i]=-99;     py_prop_debug_[i]=-99;     pz_prop_debug_[i]=-99;
    epxpx_prop_debug_[i]=-99;  epypy_prop_debug_[i]=-99;  epzpz_prop_debug_[i]=-99;
    pt_prop_debug_[i]=-99;     phi_prop_debug_[i]=-99;    eta_prop_debug_[i]=-99;
    ept_prop_debug_[i]=-99;    ephi_prop_debug_[i]=-99;   eeta_prop_debug_[i]=-99;
    invpt_prop_debug_[i] =-99; theta_prop_debug_[i] =-99;
    einvpt_prop_debug_[i]=-99; etheta_prop_debug_[i]=-99;

    // reset chi2
    layer_chi2_debug_[i]=-99;
    chi2_debug_[i]=-99;

    // reset update info
    layer_up_debug_[i]=99;
    x_up_debug_[i]=-99;      y_up_debug_[i]=-99;      z_up_debug_[i]=-99;
    exx_up_debug_[i]=-99;    eyy_up_debug_[i]=-99;    ezz_up_debug_[i]=-99;
    px_up_debug_[i]=-99;     py_up_debug_[i]=-99;     pz_up_debug_[i]=-99;
    epxpx_up_debug_[i]=-99;  epypy_up_debug_[i]=-99;  epzpz_up_debug_[i]=-99;
    pt_up_debug_[i]=-99;     phi_up_debug_[i]=-99;    eta_up_debug_[i]=-99;
    ept_up_debug_[i]=-99;    ephi_up_debug_[i]=-99;   eeta_up_debug_[i]=-99;
    invpt_up_debug_[i] =-99; theta_up_debug_[i] =-99;
    einvpt_up_debug_[i]=-99; etheta_up_debug_[i]=-99;

    // reset eta/phi bin info
    ebhit_debug_[i] = -99; ebp_debug_[i] = -99; ebm_debug_[i] = -99;
    pbhit_debug_[i] = -99; pbp_debug_[i] = -99; pbm_debug_[i] = -99;
  }
}

void TTreeValidation::fillDebugTree(const Event& ev){
  std::lock_guard<std::mutex> locker(glock_);
  resetDebugTreeArrays();
  
  auto& simtrack   = ev.simTracks_[0];       // since it is the only one!
  auto& seedtrack  = ev.seedTracks_[0];      // since it is the only one! 
  auto& buildtrack = ev.candidateTracks_[0]; // since it is the only one! 

  int mcID   = simtrack.label();
  int seedID = ev.seedTracksExtra_[seedtrack.label()].seedID(); 

  event_debug_ = ev.evtID();
  nHits_debug_ = buildtrack.nFoundHits();

  recocharge_debug_ = seedtrack.charge(); 
  mccharge_debug_   = simtrack.charge();

  // Set MC First
  pt_gen_debug_  = simtrack.pT(); 
  phi_gen_debug_ = simtrack.momPhi(); 
  eta_gen_debug_ = simtrack.momEta(); 

  x_gen_debug_   = simtrack.x(); 
  y_gen_debug_   = simtrack.y(); 
  z_gen_debug_   = simtrack.z(); 

  auto& simhits = simtrack.hitsVector(ev.layerHits_);
  for (int i = 0; i < simhits.size(); i++){ // assume one hit for layer for sim tracks...
    layer_mc_debug_[i] = i;

    x_hit_debug_[i]   = simhits[i].x();
    y_hit_debug_[i]   = simhits[i].y();
    z_hit_debug_[i]   = simhits[i].z();
    exx_hit_debug_[i] = simhits[i].exx();
    eyy_hit_debug_[i] = simhits[i].eyy();
    ezz_hit_debug_[i] = simhits[i].ezz();

    // which eta/phi bin the hit belongs to
    ebhit_debug_[i] = getEtaPartition(simhits[i].eta());
    pbhit_debug_[i] = getPhiPartition(simhits[i].phi());

    const TrackState & mcstate = simTkTSVecMap_[mcID][i];
    x_mc_debug_[i]   = mcstate.x();
    y_mc_debug_[i]   = mcstate.y();
    z_mc_debug_[i]   = mcstate.z();
    pt_mc_debug_[i]  = mcstate.pT();
    phi_mc_debug_[i] = mcstate.momPhi();
    eta_mc_debug_[i] = mcstate.momEta();

    px_mc_debug_[i]  = mcstate.px();
    py_mc_debug_[i]  = mcstate.py();
    pz_mc_debug_[i]  = mcstate.pz();
      
    invpt_mc_debug_[i] = mcstate.invpT();
    theta_mc_debug_[i] = mcstate.theta();
  }

  // CF next
  const TrackState & cfSeedTS = seedTkCFMap_[seedID];
  x_cf_debug_   = cfSeedTS.x();
  y_cf_debug_   = cfSeedTS.y();
  z_cf_debug_   = cfSeedTS.z();
  exx_cf_debug_ = cfSeedTS.exx();
  eyy_cf_debug_ = cfSeedTS.eyy();
  ezz_cf_debug_ = cfSeedTS.ezz();

  px_cf_debug_    = cfSeedTS.px();
  py_cf_debug_    = cfSeedTS.py();
  pz_cf_debug_    = cfSeedTS.pz();
  epxpx_cf_debug_ = cfSeedTS.epxpx();
  epypy_cf_debug_ = cfSeedTS.epypy();
  epzpz_cf_debug_ = cfSeedTS.epzpz();

  pt_cf_debug_   = cfSeedTS.pT();
  phi_cf_debug_  = cfSeedTS.momPhi();
  eta_cf_debug_  = cfSeedTS.momEta();
  ept_cf_debug_  = cfSeedTS.epT();
  ephi_cf_debug_ = cfSeedTS.emomPhi();
  eeta_cf_debug_ = cfSeedTS.emomEta();

  invpt_cf_debug_  = cfSeedTS.invpT();
  einvpt_cf_debug_ = cfSeedTS.einvpT();
  theta_cf_debug_  = cfSeedTS.theta();
  etheta_cf_debug_ = cfSeedTS.etheta();

  // prop states
  for (int i = 0; i < propTSLayerPairVec_.size(); i++){
    int layer = propTSLayerPairVec_[i].first; 
    layer_prop_debug_[layer] = layer;
	
    const TrackState & propTS = propTSLayerPairVec_[i].second;
    x_prop_debug_[layer]   = propTS.x();
    y_prop_debug_[layer]   = propTS.y();
    z_prop_debug_[layer]   = propTS.z();
    exx_prop_debug_[layer] = propTS.exx();
    eyy_prop_debug_[layer] = propTS.eyy();
    ezz_prop_debug_[layer] = propTS.ezz();

    px_prop_debug_[layer]    = propTS.px();
    py_prop_debug_[layer]    = propTS.py();
    pz_prop_debug_[layer]    = propTS.pz();
    epxpx_prop_debug_[layer] = propTS.epxpx();
    epypy_prop_debug_[layer] = propTS.epypy();
    epzpz_prop_debug_[layer] = propTS.epzpz();

    pt_prop_debug_[layer]   = propTS.pT();
    phi_prop_debug_[layer]  = propTS.momPhi();
    eta_prop_debug_[layer]  = propTS.momEta();
    ept_prop_debug_[layer]  = propTS.epT();
    ephi_prop_debug_[layer] = propTS.emomPhi();
    eeta_prop_debug_[layer] = propTS.emomEta();

    invpt_prop_debug_[layer]  = propTS.invpT();
    einvpt_prop_debug_[layer] = propTS.einvpT();
    theta_prop_debug_[layer]  = propTS.theta();
    etheta_prop_debug_[layer] = propTS.etheta();
  }
    
  // do the chi2
  for (int i = 0; i < chi2LayerPairVec_.size(); i++){
    int layer = chi2LayerPairVec_[i].first; 
    layer_chi2_debug_[layer] = layer;
    chi2_debug_[layer]       = chi2LayerPairVec_[i].second; 
  }

  // updated states 
  for (int i = 0; i < upTSLayerPairVec_.size(); i++){
    int layer = upTSLayerPairVec_[i].first; 
    layer_up_debug_[layer] = layer;
	
    const TrackState & upTS = upTSLayerPairVec_[i].second;
    x_up_debug_[layer]   = upTS.x();
    y_up_debug_[layer]   = upTS.y();
    z_up_debug_[layer]   = upTS.z();
    exx_up_debug_[layer] = upTS.exx();
    eyy_up_debug_[layer] = upTS.eyy();
    ezz_up_debug_[layer] = upTS.ezz();

    px_up_debug_[layer]    = upTS.px();
    py_up_debug_[layer]    = upTS.py();
    pz_up_debug_[layer]    = upTS.pz();
    epxpx_up_debug_[layer] = upTS.epxpx();
    epypy_up_debug_[layer] = upTS.epypy();
    epzpz_up_debug_[layer] = upTS.epzpz();

    pt_up_debug_[layer]   = upTS.pT();
    phi_up_debug_[layer]  = upTS.momPhi();
    eta_up_debug_[layer]  = upTS.momEta();
    ept_up_debug_[layer]  = upTS.epT();
    ephi_up_debug_[layer] = upTS.emomPhi();
    eeta_up_debug_[layer] = upTS.emomEta();

    invpt_up_debug_[layer]  = upTS.invpT();
    einvpt_up_debug_[layer] = upTS.einvpT();
    theta_up_debug_[layer]  = upTS.theta();
    etheta_up_debug_[layer] = upTS.etheta();
  }

  // see what it predicted! (reuse branchval)
  for (TkIDToBVVMMIter seediter = seedToBranchValVecLayMapMap_.begin(); seediter != seedToBranchValVecLayMapMap_.end(); ++seediter){
    for (BVVLMiter layiter = (*seediter).second.begin(); layiter != (*seediter).second.end(); ++layiter){
      const auto& BranchValVec((*layiter).second);
      const int cands = BranchValVec.size();
      int layer = (*layiter).first; // first index here is layer
      for (int cand = 0; cand < cands; cand++){ // loop over input candidates at this layer for this seed
	const auto& BranchVal(BranchValVec[cand]); // grab the branch validation object
	
	ebp_debug_[layer]  = BranchVal.etaBinPlus;
	ebm_debug_[layer]  = BranchVal.etaBinMinus;
	pbp_debug_[layer]  = BranchVal.phiBinPlus;
	pbm_debug_[layer]  = BranchVal.phiBinMinus;

      }
    }
  }

  // fill it once per track (i.e. once per track by definition)
  debugtree_->Fill();
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

void TTreeValidation::fillSeedInfoTree(const TripletIdxVec& hit_triplets, const Event& ev) {
  evtID_seedinfo_ = ev.evtID();

  const auto & evt_lay_hits = ev.layerHits_;

  for (auto&& hitTriplet : hit_triplets) {
    mcID_seedinfo_    = ev.simHitsInfo_[evt_lay_hits[0][hitTriplet[0]].mcHitID()].mcTrackID();
    if (  ev.simHitsInfo_[evt_lay_hits[1][hitTriplet[1]].mcHitID()].mcTrackID() == mcID_seedinfo_){
      if (ev.simHitsInfo_[evt_lay_hits[2][hitTriplet[2]].mcHitID()].mcTrackID() == mcID_seedinfo_){
	pt_gen_seedinfo_  = ev.simTracks_[mcID_seedinfo_].pT();
	eta_gen_seedinfo_ = ev.simTracks_[mcID_seedinfo_].momEta();
	phi_gen_seedinfo_ = ev.simTracks_[mcID_seedinfo_].momPhi();

	const float x0 = evt_lay_hits[0][hitTriplet[0]].x();
	const float y0 = evt_lay_hits[0][hitTriplet[0]].y();
	const float x1 = evt_lay_hits[1][hitTriplet[1]].x();
	const float y1 = evt_lay_hits[1][hitTriplet[1]].y();
	const float x2 = evt_lay_hits[2][hitTriplet[2]].x();
	const float y2 = evt_lay_hits[2][hitTriplet[2]].y();

	// now fit a circle, extract pT and d0 from center and radius
	const float mr = (y1-y0)/(x1-x0);
	const float mt = (y2-y1)/(x2-x1);

	a  = (mr*mt*(y2-y0) + mr*(x1+x2) - mt*(x0+x1))/(2.*(mr-mt));
        b  = -1.*(a-(x0+x1)/2.)/mr + (y0+y1)/2.;
	r  = getHypot(x0-a,y0-b);
	d0 = getHypot(a,b)-r;
	pass = false;
	if ((r >= Config::maxCurvR) && (fabs(d0) <= 0.1)) {
	  pass = true;
	} // d0 cut 1mm, pT cut 0.5 GeV

	seedinfotree_->Fill();
      }
    }
  }
}

void TTreeValidation::fillSeedTree(const TripletIdxVec& hit_triplets, const TripletIdxVec& filtered_triplets, const Event& ev) {
  evtID_seed_ = ev.evtID();
  const auto & evt_lay_hits = ev.layerHits_;

  int correct_all = 0;
  for (auto&& hitTriplet : hit_triplets){
    if (  ev.simHitsInfo_[evt_lay_hits[0][hitTriplet[0]].mcHitID()].mcTrackID() == ev.simHitsInfo_[evt_lay_hits[1][hitTriplet[1]].mcHitID()].mcTrackID()){
      if (ev.simHitsInfo_[evt_lay_hits[2][hitTriplet[2]].mcHitID()].mcTrackID() == ev.simHitsInfo_[evt_lay_hits[1][hitTriplet[1]].mcHitID()].mcTrackID()){
	correct_all++;
      }
    }
  }
  nTkAll_   = hit_triplets.size();
  nTkAllMC_ = correct_all; 

  int correct_cut = 0;
  for (auto&& hitTriplet : filtered_triplets){
    if (  ev.simHitsInfo_[evt_lay_hits[0][hitTriplet[0]].mcHitID()].mcTrackID() == ev.simHitsInfo_[evt_lay_hits[1][hitTriplet[1]].mcHitID()].mcTrackID()){
      if (ev.simHitsInfo_[evt_lay_hits[2][hitTriplet[2]].mcHitID()].mcTrackID() == ev.simHitsInfo_[evt_lay_hits[1][hitTriplet[1]].mcHitID()].mcTrackID()){
	correct_cut++;
      }
    }
  }
  nTkCut_   = filtered_triplets.size();
  nTkCutMC_ = correct_cut; 

  seedtree_->Fill();
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
      const int layer = seedtrack.foundLayers().back(); // last layer seed ended up on
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

      hitchi2_seed_eff_   = seedtrack.chi2(); // currently not being used
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
      const int layer = buildtrack.foundLayers().back(); // last layer build ended up on
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
      const int layer = fittrack.foundLayers().front(); // last layer fit ended up on
#else
      const int layer = fittrack.foundLayers().back(); // last layer fit ended up on
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

    hitchi2_seed_FR_ = seedtrack.chi2(); //--> not currently used

    // sim info for seed track
    mcID_seed_FR_ = seedextra.mcTrackID();
    if (mcID_seed_FR_ != 999999){ // store sim info at that final layer!!! --> gen info stored only in eff tree
      mcmask_seed_FR_ = 1; // matched track to sim

      const int layer = seedtrack.foundLayers().back(); // last layer fit ended up on
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
	
	const int layer = buildtrack.foundLayers().back(); // last layer fit ended up on
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
	const int layer = fittrack.foundLayers().front(); // last layer fiting ended up on
#else
	const int layer = fittrack.foundLayers().back(); // last layer fiting ended up on
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
      const int cflayer = seedtrack.foundLayers().front(); //  layer for which cf parameters are calculated with respect to (either first or last layer of track!)
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
      const int cflayer = fittrack.foundLayers().back(); // last layer fit ended up on
#else
      const int cflayer = fittrack.foundLayers().front(); // last layer fit ended up on
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
