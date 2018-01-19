#include "TTreeValidation.h"
#include "Event.h"
#include "Config.h"
#include "Propagation.h"
#ifndef NO_ROOT

TTreeValidation::TTreeValidation(std::string fileName)
{
  std::lock_guard<std::mutex> locker(glock_);
  gROOT->ProcessLine("#include <vector>");
  f_ = TFile::Open(fileName.c_str(), "recreate");

  if (Config::root_val) 
  { 
    TTreeValidation::initializeEfficiencyTree();
    TTreeValidation::initializeFakeRateTree();  
  }
  if (Config::cmssw_val)
  {
    TTreeValidation::initializeCMSSWEfficiencyTree();
    TTreeValidation::initializeCMSSWFakeRateTree();
  }
  if (Config::fit_val) 
  {
    for (int i = 0; i < nfvs_; ++i) fvs_[i].resize(Config::nTotalLayers);
    TTreeValidation::initializeFitTree();
  }
  TTreeValidation::initializeConfigTree();
}

TTreeValidation::~TTreeValidation()
{
  if (Config::root_val) 
  {
    delete efftree_;
    delete frtree_;
  }
  if (Config::cmssw_val)
  {
    delete cmsswefftree_;
    delete cmsswfrtree_;
  }
  if (Config::fit_val) 
  {
    delete fittree_;
  }
  delete configtree_;
  delete f_;
}

void TTreeValidation::initializeEfficiencyTree()
{  
  // efficiency validation
  efftree_ = new TTree("efftree","efftree");
  efftree_->Branch("evtID",&evtID_eff_);
  efftree_->Branch("mcID",&mcID_eff_);

  efftree_->Branch("nHits_mc",&nHits_mc_eff_);
  efftree_->Branch("lastlyr_mc",&lastlyr_mc_eff_);

  efftree_->Branch("seedID_seed",&seedID_seed_eff_);
  efftree_->Branch("seedID_build",&seedID_build_eff_);
  efftree_->Branch("seedID_fit",&seedID_fit_eff_);

  efftree_->Branch("x_mc_gen",&x_mc_gen_eff_);
  efftree_->Branch("y_mc_gen",&y_mc_gen_eff_);
  efftree_->Branch("z_mc_gen",&z_mc_gen_eff_);

  efftree_->Branch("pt_mc_gen",&pt_mc_gen_eff_);
  efftree_->Branch("phi_mc_gen",&phi_mc_gen_eff_);
  efftree_->Branch("eta_mc_gen",&eta_mc_gen_eff_);

  efftree_->Branch("mcmask_seed",&mcmask_seed_eff_);
  efftree_->Branch("mcmask_build",&mcmask_build_eff_);
  efftree_->Branch("mcmask_fit",&mcmask_fit_eff_);

  efftree_->Branch("xhit_seed",&xhit_seed_eff_);
  efftree_->Branch("xhit_build",&xhit_build_eff_);
  efftree_->Branch("xhit_fit",&xhit_fit_eff_);

  efftree_->Branch("yhit_seed",&yhit_seed_eff_);
  efftree_->Branch("yhit_build",&yhit_build_eff_);
  efftree_->Branch("yhit_fit",&yhit_fit_eff_);

  efftree_->Branch("zhit_seed",&zhit_seed_eff_);
  efftree_->Branch("zhit_build",&zhit_build_eff_);
  efftree_->Branch("zhit_fit",&zhit_fit_eff_);

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

  efftree_->Branch("lastlyr_seed",&lastlyr_seed_eff_);
  efftree_->Branch("lastlyr_build",&lastlyr_build_eff_);
  efftree_->Branch("lastlyr_fit",&lastlyr_fit_eff_);

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

void TTreeValidation::initializeFakeRateTree()
{
  // fake rate validation
  frtree_ = new TTree("frtree","frtree");

  frtree_->Branch("evtID",&evtID_FR_);
  frtree_->Branch("seedID",&seedID_FR_);

  frtree_->Branch("seedmask_seed",&seedmask_seed_FR_);
  frtree_->Branch("seedmask_build",&seedmask_build_FR_);
  frtree_->Branch("seedmask_fit",&seedmask_fit_FR_);

  frtree_->Branch("xhit_seed",&xhit_seed_FR_);
  frtree_->Branch("xhit_build",&xhit_build_FR_);
  frtree_->Branch("xhit_fit",&xhit_fit_FR_);

  frtree_->Branch("yhit_seed",&yhit_seed_FR_);
  frtree_->Branch("yhit_build",&yhit_build_FR_);
  frtree_->Branch("yhit_fit",&yhit_fit_FR_);

  frtree_->Branch("zhit_seed",&zhit_seed_FR_);
  frtree_->Branch("zhit_build",&zhit_build_FR_);
  frtree_->Branch("zhit_fit",&zhit_fit_FR_);

  frtree_->Branch("pt_seed",&pt_seed_FR_);
  frtree_->Branch("ept_seed",&ept_seed_FR_);
  frtree_->Branch("pt_build",&pt_build_FR_);
  frtree_->Branch("ept_build",&ept_build_FR_);
  frtree_->Branch("pt_fit",&pt_fit_FR_);
  frtree_->Branch("ept_fit",&ept_fit_FR_);

  frtree_->Branch("phi_seed",&phi_seed_FR_);
  frtree_->Branch("ephi_seed",&ephi_seed_FR_);
  frtree_->Branch("phi_build",&phi_build_FR_);
  frtree_->Branch("ephi_build",&ephi_build_FR_);
  frtree_->Branch("phi_fit",&phi_fit_FR_);
  frtree_->Branch("ephi_fit",&ephi_fit_FR_);

  frtree_->Branch("eta_seed",&eta_seed_FR_);
  frtree_->Branch("eeta_seed",&eeta_seed_FR_);
  frtree_->Branch("eta_build",&eta_build_FR_);
  frtree_->Branch("eeta_build",&eeta_build_FR_);
  frtree_->Branch("eta_fit",&eta_fit_FR_);
  frtree_->Branch("eeta_fit",&eeta_fit_FR_);

  frtree_->Branch("nHits_seed",&nHits_seed_FR_);
  frtree_->Branch("nHits_build",&nHits_build_FR_);
  frtree_->Branch("nHits_fit",&nHits_fit_FR_);

  frtree_->Branch("nHitsMatched_seed",&nHitsMatched_seed_FR_);
  frtree_->Branch("nHitsMatched_build",&nHitsMatched_build_FR_);
  frtree_->Branch("nHitsMatched_fit",&nHitsMatched_fit_FR_);

  frtree_->Branch("fracHitsMatched_seed",&fracHitsMatched_seed_FR_);
  frtree_->Branch("fracHitsMatched_build",&fracHitsMatched_build_FR_);
  frtree_->Branch("fracHitsMatched_fit",&fracHitsMatched_fit_FR_);

  frtree_->Branch("lastlyr_seed",&lastlyr_seed_FR_);
  frtree_->Branch("lastlyr_build",&lastlyr_build_FR_);
  frtree_->Branch("lastlyr_fit",&lastlyr_fit_FR_);

  frtree_->Branch("hitchi2_seed",&hitchi2_seed_FR_);
  frtree_->Branch("hitchi2_build",&hitchi2_build_FR_);
  frtree_->Branch("hitchi2_fit",&hitchi2_fit_FR_);

  // sim info of seed,build,fit tracks
  frtree_->Branch("mcID_seed",&mcID_seed_FR_);
  frtree_->Branch("mcID_build",&mcID_build_FR_);
  frtree_->Branch("mcID_fit",&mcID_fit_FR_);
  
  frtree_->Branch("mcmask_seed",&mcmask_seed_FR_);
  frtree_->Branch("mcmask_build",&mcmask_build_FR_);
  frtree_->Branch("mcmask_fit",&mcmask_fit_FR_);

  frtree_->Branch("pt_mc_seed",&pt_mc_seed_FR_);
  frtree_->Branch("pt_mc_build",&pt_mc_build_FR_);
  frtree_->Branch("pt_mc_fit",&pt_mc_fit_FR_);

  frtree_->Branch("phi_mc_seed",&phi_mc_seed_FR_);
  frtree_->Branch("phi_mc_build",&phi_mc_build_FR_);
  frtree_->Branch("phi_mc_fit",&phi_mc_fit_FR_);

  frtree_->Branch("eta_mc_seed",&eta_mc_seed_FR_);
  frtree_->Branch("eta_mc_build",&eta_mc_build_FR_);
  frtree_->Branch("eta_mc_fit",&eta_mc_fit_FR_);

  frtree_->Branch("nHits_mc_seed",&nHits_mc_seed_FR_);
  frtree_->Branch("nHits_mc_build",&nHits_mc_build_FR_);
  frtree_->Branch("nHits_mc_fit",&nHits_mc_fit_FR_);

  frtree_->Branch("lastlyr_mc_seed",&lastlyr_mc_seed_FR_);
  frtree_->Branch("lastlyr_mc_build",&lastlyr_mc_build_FR_);
  frtree_->Branch("lastlyr_mc_fit",&lastlyr_mc_fit_FR_);

  frtree_->Branch("helixchi2_seed",&helixchi2_seed_FR_);
  frtree_->Branch("helixchi2_build",&helixchi2_build_FR_);
  frtree_->Branch("helixchi2_fit",&helixchi2_fit_FR_);

  frtree_->Branch("duplmask_seed",&duplmask_seed_FR_);
  frtree_->Branch("duplmask_build",&duplmask_build_FR_);
  frtree_->Branch("duplmask_fit",&duplmask_fit_FR_);

  frtree_->Branch("iTkMatches_seed",&iTkMatches_seed_FR_);
  frtree_->Branch("iTkMatches_build",&iTkMatches_build_FR_);
  frtree_->Branch("iTkMatches_fit",&iTkMatches_fit_FR_);
}

void TTreeValidation::initializeConfigTree()
{
  // include config ++ real seeding parameters ...
  configtree_ = new TTree("configtree","configtree");

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

void TTreeValidation::initializeCMSSWEfficiencyTree()
{  
  // cmssw reco track efficiency validation
  cmsswefftree_ = new TTree("cmsswefftree","cmsswefftree");
  cmsswefftree_->Branch("evtID",&evtID_ceff_);
  cmsswefftree_->Branch("cmsswID",&cmsswID_ceff_);
  cmsswefftree_->Branch("seedID_cmssw",&seedID_cmssw_ceff_);

  // CMSSW
  cmsswefftree_->Branch("x_cmssw",&x_cmssw_ceff_);
  cmsswefftree_->Branch("y_cmssw",&y_cmssw_ceff_);
  cmsswefftree_->Branch("z_cmssw",&z_cmssw_ceff_);

  cmsswefftree_->Branch("pt_cmssw",&pt_cmssw_ceff_);
  cmsswefftree_->Branch("phi_cmssw",&phi_cmssw_ceff_);
  cmsswefftree_->Branch("eta_cmssw",&eta_cmssw_ceff_);

  cmsswefftree_->Branch("nHits_cmssw",&nHits_cmssw_ceff_);
  cmsswefftree_->Branch("nLayers_cmssw",&nLayers_cmssw_ceff_);
  cmsswefftree_->Branch("lastlyr_cmssw",&lastlyr_cmssw_ceff_);

  // Build
  cmsswefftree_->Branch("cmsswmask_build",&cmsswmask_build_ceff_);
  cmsswefftree_->Branch("seedID_build",&seedID_build_ceff_);
  cmsswefftree_->Branch("mcTrackID_build",&mcTrackID_build_ceff_);

  cmsswefftree_->Branch("pt_build",&pt_build_ceff_);
  cmsswefftree_->Branch("ept_build",&ept_build_ceff_);
  cmsswefftree_->Branch("phi_build",&phi_build_ceff_);
  cmsswefftree_->Branch("ephi_build",&ephi_build_ceff_);
  cmsswefftree_->Branch("eta_build",&eta_build_ceff_);
  cmsswefftree_->Branch("eeta_build",&eeta_build_ceff_);

  cmsswefftree_->Branch("nHits_build",&nHits_build_ceff_);
  cmsswefftree_->Branch("nLayers_build",&nLayers_build_ceff_);
  cmsswefftree_->Branch("nHitsMatched_build",&nHitsMatched_build_ceff_);
  cmsswefftree_->Branch("fracHitsMatched_build",&fracHitsMatched_build_ceff_);
  cmsswefftree_->Branch("lastlyr_build",&lastlyr_build_ceff_);

  cmsswefftree_->Branch("xhit_build",&xhit_build_ceff_);
  cmsswefftree_->Branch("yhit_build",&yhit_build_ceff_);
  cmsswefftree_->Branch("zhit_build",&zhit_build_ceff_);

  cmsswefftree_->Branch("hitchi2_build",&hitchi2_build_ceff_);
  cmsswefftree_->Branch("helixchi2_build",&helixchi2_build_ceff_);
  cmsswefftree_->Branch("dphi_build",&dphi_build_ceff_);

  cmsswefftree_->Branch("duplmask_build",&duplmask_build_ceff_);
  cmsswefftree_->Branch("nTkMatches_build",&nTkMatches_build_ceff_);

  // Fit
  cmsswefftree_->Branch("cmsswmask_fit",&cmsswmask_fit_ceff_);
  cmsswefftree_->Branch("seedID_fit",&seedID_fit_ceff_);
  cmsswefftree_->Branch("mcTrackID_fit",&mcTrackID_fit_ceff_);

  cmsswefftree_->Branch("pt_fit",&pt_fit_ceff_);
  cmsswefftree_->Branch("ept_fit",&ept_fit_ceff_);
  cmsswefftree_->Branch("phi_fit",&phi_fit_ceff_);
  cmsswefftree_->Branch("ephi_fit",&ephi_fit_ceff_);
  cmsswefftree_->Branch("eta_fit",&eta_fit_ceff_);
  cmsswefftree_->Branch("eeta_fit",&eeta_fit_ceff_);

  cmsswefftree_->Branch("nHits_fit",&nHits_fit_ceff_);
  cmsswefftree_->Branch("nLayers_fit",&nLayers_fit_ceff_);
  cmsswefftree_->Branch("nHitsMatched_fit",&nHitsMatched_fit_ceff_);
  cmsswefftree_->Branch("fracHitsMatched_fit",&fracHitsMatched_fit_ceff_);
  cmsswefftree_->Branch("lastlyr_fit",&lastlyr_fit_ceff_);

  cmsswefftree_->Branch("xhit_fit",&xhit_fit_ceff_);
  cmsswefftree_->Branch("yhit_fit",&yhit_fit_ceff_);
  cmsswefftree_->Branch("zhit_fit",&zhit_fit_ceff_);

  cmsswefftree_->Branch("hitchi2_fit",&hitchi2_fit_ceff_);
  cmsswefftree_->Branch("helixchi2_fit",&helixchi2_fit_ceff_);
  cmsswefftree_->Branch("dphi_fit",&dphi_fit_ceff_);

  cmsswefftree_->Branch("duplmask_fit",&duplmask_fit_ceff_);
  cmsswefftree_->Branch("nTkMatches_fit",&nTkMatches_fit_ceff_);
}

void TTreeValidation::initializeCMSSWFakeRateTree()
{  
  // cmssw reco track efficiency validation
  cmsswfrtree_ = new TTree("cmsswfrtree","cmsswfrtree");

  cmsswfrtree_->Branch("evtID",&evtID_cFR_);
  cmsswfrtree_->Branch("seedID",&seedID_cFR_);
  cmsswfrtree_->Branch("mcTrackID",&mcTrackID_cFR_);

  // build
  cmsswfrtree_->Branch("cmsswID_build",&cmsswID_build_cFR_);
  cmsswfrtree_->Branch("cmsswmask_build",&cmsswmask_build_cFR_);

  cmsswfrtree_->Branch("pt_build",&pt_build_cFR_);
  cmsswfrtree_->Branch("ept_build",&ept_build_cFR_);
  cmsswfrtree_->Branch("phi_build",&phi_build_cFR_);
  cmsswfrtree_->Branch("ephi_build",&ephi_build_cFR_);
  cmsswfrtree_->Branch("eta_build",&eta_build_cFR_);
  cmsswfrtree_->Branch("eeta_build",&eeta_build_cFR_);

  cmsswfrtree_->Branch("nHits_build",&nHits_build_cFR_);
  cmsswfrtree_->Branch("nLayers_build",&nLayers_build_cFR_);
  cmsswfrtree_->Branch("nHitsMatched_build",&nHitsMatched_build_cFR_);
  cmsswfrtree_->Branch("fracHitsMatched_build",&fracHitsMatched_build_cFR_);
  cmsswfrtree_->Branch("lastlyr_build",&lastlyr_build_cFR_);

  cmsswfrtree_->Branch("xhit_build",&xhit_build_cFR_);
  cmsswfrtree_->Branch("yhit_build",&yhit_build_cFR_);
  cmsswfrtree_->Branch("zhit_build",&zhit_build_cFR_);

  cmsswfrtree_->Branch("hitchi2_build",&hitchi2_build_cFR_);
  cmsswfrtree_->Branch("helixchi2_build",&helixchi2_build_cFR_);
  cmsswfrtree_->Branch("dphi_build",&dphi_build_cFR_);

  cmsswfrtree_->Branch("duplmask_build",&duplmask_build_cFR_);
  cmsswfrtree_->Branch("iTkMatches_build",&iTkMatches_build_cFR_);

  cmsswfrtree_->Branch("seedID_cmssw_build",&seedID_cmssw_build_cFR_);

  cmsswfrtree_->Branch("x_cmssw_build",&x_cmssw_build_cFR_);
  cmsswfrtree_->Branch("y_cmssw_build",&y_cmssw_build_cFR_);
  cmsswfrtree_->Branch("z_cmssw_build",&z_cmssw_build_cFR_);

  cmsswfrtree_->Branch("pt_cmssw_build",&pt_cmssw_build_cFR_);
  cmsswfrtree_->Branch("phi_cmssw_build",&phi_cmssw_build_cFR_);
  cmsswfrtree_->Branch("eta_cmssw_build",&eta_cmssw_build_cFR_);

  cmsswfrtree_->Branch("nHits_cmssw_build",&nHits_cmssw_build_cFR_);
  cmsswfrtree_->Branch("nLayers_cmssw_build",&nLayers_cmssw_build_cFR_);
  cmsswfrtree_->Branch("lastlyr_cmssw_build",&lastlyr_cmssw_build_cFR_);

  // fit
  cmsswfrtree_->Branch("cmsswID_fit",&cmsswID_fit_cFR_);
  cmsswfrtree_->Branch("cmsswmask_fit",&cmsswmask_fit_cFR_);

  cmsswfrtree_->Branch("pt_fit",&pt_fit_cFR_);
  cmsswfrtree_->Branch("ept_fit",&ept_fit_cFR_);
  cmsswfrtree_->Branch("phi_fit",&phi_fit_cFR_);
  cmsswfrtree_->Branch("ephi_fit",&ephi_fit_cFR_);
  cmsswfrtree_->Branch("eta_fit",&eta_fit_cFR_);
  cmsswfrtree_->Branch("eeta_fit",&eeta_fit_cFR_);

  cmsswfrtree_->Branch("nHits_fit",&nHits_fit_cFR_);
  cmsswfrtree_->Branch("nLayers_fit",&nLayers_fit_cFR_);
  cmsswfrtree_->Branch("nHitsMatched_fit",&nHitsMatched_fit_cFR_);
  cmsswfrtree_->Branch("fracHitsMatched_fit",&fracHitsMatched_fit_cFR_);
  cmsswfrtree_->Branch("lastlyr_fit",&lastlyr_fit_cFR_);

  cmsswfrtree_->Branch("xhit_fit",&xhit_fit_cFR_);
  cmsswfrtree_->Branch("yhit_fit",&yhit_fit_cFR_);
  cmsswfrtree_->Branch("zhit_fit",&zhit_fit_cFR_);

  cmsswfrtree_->Branch("hitchi2_fit",&hitchi2_fit_cFR_);
  cmsswfrtree_->Branch("helixchi2_fit",&helixchi2_fit_cFR_);
  cmsswfrtree_->Branch("dphi_fit",&dphi_fit_cFR_);

  cmsswfrtree_->Branch("duplmask_fit",&duplmask_fit_cFR_);
  cmsswfrtree_->Branch("iTkMatches_fit",&iTkMatches_fit_cFR_);

  cmsswfrtree_->Branch("seedID_cmssw_fit",&seedID_cmssw_fit_cFR_);

  cmsswfrtree_->Branch("x_cmssw_fit",&x_cmssw_fit_cFR_);
  cmsswfrtree_->Branch("y_cmssw_fit",&y_cmssw_fit_cFR_);
  cmsswfrtree_->Branch("z_cmssw_fit",&z_cmssw_fit_cFR_);

  cmsswfrtree_->Branch("pt_cmssw_fit",&pt_cmssw_fit_cFR_);
  cmsswfrtree_->Branch("phi_cmssw_fit",&phi_cmssw_fit_cFR_);
  cmsswfrtree_->Branch("eta_cmssw_fit",&eta_cmssw_fit_cFR_);

  cmsswfrtree_->Branch("nHits_cmssw_fit",&nHits_cmssw_fit_cFR_);
  cmsswfrtree_->Branch("nLayers_cmssw_fit",&nLayers_cmssw_fit_cFR_);
  cmsswfrtree_->Branch("lastlyr_cmssw_fit",&lastlyr_cmssw_fit_cFR_);
}

void TTreeValidation::initializeFitTree()
{
  ntotallayers_fit_ = Config::nTotalLayers;
  
  fittree_ = new TTree("fittree","fittree");

  fittree_->Branch("ntotallayers",&ntotallayers_fit_,"ntotallayers_fit_/I");
  fittree_->Branch("tkid",&tkid_fit_,"tkid_fit_/I");
  fittree_->Branch("evtid",&evtid_fit_,"evtid_fit_/I");

  fittree_->Branch("z_prop",&z_prop_fit_,"z_prop_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("ez_prop",&ez_prop_fit_,"ez_prop_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("z_hit",&z_hit_fit_,"z_hit_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("ez_hit",&ez_hit_fit_,"ez_hit_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("z_sim",&z_sim_fit_,"z_sim_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("ez_sim",&ez_sim_fit_,"ez_sim_fit_[ntotallayers_fit_]/F");

  fittree_->Branch("pphi_prop",&pphi_prop_fit_,"pphi_prop_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("epphi_prop",&epphi_prop_fit_,"epphi_prop_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("pphi_hit",&pphi_hit_fit_,"pphi_hit_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("epphi_hit",&epphi_hit_fit_,"epphi_hit_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("pphi_sim",&pphi_sim_fit_,"pphi_sim_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("epphi_sim",&epphi_sim_fit_,"epphi_sim_fit_[ntotallayers_fit_]/F");

  fittree_->Branch("pt_up",&pt_up_fit_,"pt_up_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("ept_up",&ept_up_fit_,"ept_up_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("pt_sim",&pt_sim_fit_,"pt_sim_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("ept_sim",&ept_sim_fit_,"ept_sim_fit_[ntotallayers_fit_]/F");

  fittree_->Branch("mphi_up",&mphi_up_fit_,"mphi_up_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("emphi_up",&emphi_up_fit_,"emphi_up_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("mphi_sim",&mphi_sim_fit_,"mphi_sim_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("emphi_sim",&emphi_sim_fit_,"emphi_sim_fit_[ntotallayers_fit_]/F");

  fittree_->Branch("meta_up",&meta_up_fit_,"meta_up_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("emeta_up",&emeta_up_fit_,"emeta_up_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("meta_sim",&meta_sim_fit_,"meta_sim_fit_[ntotallayers_fit_]/F");
  fittree_->Branch("emeta_sim",&emeta_sim_fit_,"emeta_sim_fit_[ntotallayers_fit_]/F");
}

void TTreeValidation::alignTracks(TrackVec& evt_tracks, TrackExtraVec& evt_extras, bool alignExtra)
{
  // redo trackExtras first if necessary
  if (alignExtra)
  {
    TrackExtraVec trackExtra_tmp(evt_tracks.size());

    // align temporary tkExVec with new track collection ordering
    for (int itrack = 0; itrack < evt_tracks.size(); itrack++)
    { 
      trackExtra_tmp[itrack] = evt_extras[evt_tracks[itrack].label()]; // label is old seedID!
    }
    
    // now copy the temporary back in the old one
    evt_extras = trackExtra_tmp;
  }  

  // redo track labels to match index in vector
  for (int itrack = 0; itrack < evt_tracks.size(); itrack++)
  {
    evt_tracks[itrack].setLabel(itrack);
  }
}

void TTreeValidation::collectFitInfo(const FitVal & tmpfitval, int tkid, int layer)
{
  std::lock_guard<std::mutex> locker(glock_);

  fitValTkMapMap_[tkid][layer] = tmpfitval;
}

void TTreeValidation::resetValidationMaps()
{
  std::lock_guard<std::mutex> locker(glock_);
  // reset fit validation map
  fitValTkMapMap_.clear();
  
  // reset map of sim tracks to reco tracks
  simToSeedMap_.clear();
  simToBuildMap_.clear();
  simToFitMap_.clear();

  // reset map of seed tracks to reco tracks
  seedToBuildMap_.clear();
  seedToFitMap_.clear();

  // reset map of cmssw tracks to reco tracks
  cmsswToBuildMap_.clear();
  cmsswToFitMap_.clear();

  // reset special map of seed labels to cmssw tracks
  seedToCmsswMap_.clear();

  // reset special map of matching build tracks exactly to cmssw tracks through seedIDs
  buildToCmsswMap_.clear();

  // reset special maps used for pairing build to fit tracks CMSSW only
  buildToFitMap_.clear();
  fitToBuildMap_.clear();
}

void TTreeValidation::setTrackExtras(Event& ev)
{
  std::lock_guard<std::mutex> locker(glock_);

  const auto& layerhits = ev.layerHits_;

  if (Config::root_val)
  {
    const auto& simhits = ev.simHitsInfo_;
    const auto& simtracks = ev.simTracks_;
    const auto& seedtracks = ev.seedTracks_;
          auto& seedextras = ev.seedTracksExtra_;
    const auto& buildtracks = ev.candidateTracks_;
          auto& buildextras = ev.candidateTracksExtra_;
    const auto& fittracks = ev.fitTracks_;
          auto& fitextras = ev.fitTracksExtra_;

    // set mcTrackID for seed tracks
    for (int itrack = 0; itrack < seedtracks.size(); itrack++)
    {
      const auto& track = seedtracks[itrack];
            auto& extra = seedextras[itrack];
      extra.setMCTrackIDInfo(track, layerhits, simhits, simtracks, true); // otherwise seeds are completely unmatched in ToyMC Sim Seeds
    }
    
    // set mcTrackID for built tracks
    for (int itrack = 0; itrack < buildtracks.size(); itrack++)
    {
      const auto& track = buildtracks[itrack];
            auto& extra = buildextras[itrack];
      if (Config::seedInput == cmsswSeeds || Config::seedInput == findSeeds) {extra.setMCTrackIDInfo(track, layerhits, simhits, simtracks, false);}
      else {extra.setMCTrackIDInfoByLabel(track, layerhits, simhits, simtracks);}
    }
  
    // set mcTrackID for fit tracks
    for (int itrack = 0; itrack < fittracks.size(); itrack++)
    {
      const auto& track = fittracks[itrack];
            auto& extra = fitextras[itrack];
      if (Config::seedInput == cmsswSeeds || Config::seedInput == findSeeds) {extra.setMCTrackIDInfo(track, layerhits, simhits, simtracks, false);}
      else {extra.setMCTrackIDInfoByLabel(track, layerhits, simhits, simtracks);}
    }
  }

  if (Config::cmssw_val)
  {    
    // store mcTrackID and seedID correctly
    storeSeedAndMCID(ev);

    const auto& cmsswtracks = ev.cmsswTracks_;
          auto& cmsswextras = ev.cmsswTracksExtra_;
    const auto& buildtracks = ev.candidateTracks_;
          auto& buildextras = ev.candidateTracksExtra_;
    const auto& fittracks   = ev.fitTracks_;
          auto& fitextras   = ev.fitTracksExtra_;

    // store reduced parameters, hit map of cmssw tracks, and global hit map
    RedTrackVec reducedCMSSW;
    LayIdxIDVecMapMap cmsswHitIDMap;
    setupCMSSWMatching(ev,reducedCMSSW,cmsswHitIDMap);

    // set cmsswTrackID for built tracks
    for (int itrack = 0; itrack < buildtracks.size(); itrack++)
    {
      const auto& track = buildtracks[itrack];
            auto& extra = buildextras[itrack];

      if (Config::cmsswMatching == trkParamBased)
      {	
	extra.setCMSSWTrackIDInfoByTrkParams(track, layerhits, cmsswtracks, reducedCMSSW, false);
      }
      else if (Config::cmsswMatching == hitBased)
      {
	extra.setCMSSWTrackIDInfoByHits(track, cmsswHitIDMap, cmsswtracks, reducedCMSSW);
      }
      else if (Config::cmsswMatching == labelBased) // can only be used if using pure seeds!
      {
	extra.setCMSSWTrackIDInfoByLabel(track, layerhits, cmsswtracks, reducedCMSSW[cmsswtracks[buildToCmsswMap_[track.label()]].label()]);	
      }
      else 
      {
	std::cerr << "Specified CMSSW validation, but using an incorrect matching option! Exiting..." << std::endl;
	exit(1);
      }
    }

    // set cmsswTrackID for fit tracks --> at this point only do track param matching!
    for (int itrack = 0; itrack < fittracks.size(); itrack++)
    {
      const auto& track = fittracks[itrack];
            auto& extra = fitextras[itrack];

      extra.setCMSSWTrackIDInfoByTrkParams(track, layerhits, cmsswtracks, reducedCMSSW, true);
    }    
  }
}

void TTreeValidation::makeSimTkToRecoTksMaps(Event& ev)
{
  std::lock_guard<std::mutex> locker(glock_);
  // map sim track ids to reco labels sort by each (simTracks set in order by default!)
  TTreeValidation::mapRefTkToRecoTks(ev.seedTracks_,ev.seedTracksExtra_,simToSeedMap_);
  TTreeValidation::mapRefTkToRecoTks(ev.candidateTracks_,ev.candidateTracksExtra_,simToBuildMap_);
  TTreeValidation::mapRefTkToRecoTks(ev.fitTracks_,ev.fitTracksExtra_,simToFitMap_);
}

void TTreeValidation::mapRefTkToRecoTks(const TrackVec& evt_tracks, TrackExtraVec& evt_extras, TkIDToTkIDVecMap& refTkMap)
{
  for (auto itrack = 0; itrack < evt_tracks.size(); ++itrack)
  {
    auto&& track(evt_tracks[itrack]);
    auto&& extra(evt_extras[itrack]);
    if (Config::root_val)
    {
      if (extra.mcTrackID() >= 0) // skip fakes, don't store them at all in sim map
      {
	refTkMap[extra.mcTrackID()].push_back(track.label()); // store vector of reco tk labels, mapped to the sim track label (i.e. mcTrackID)
      }
    }
    if (Config::cmssw_val)
    {
      if (extra.cmsswTrackID() >= 0) // skip fakes, don't store them at all in cmssw map
      {
	refTkMap[extra.cmsswTrackID()].push_back(track.label()); // store vector of reco tk labels, mapped to the cmssw track label (i.e. cmsswTrackID)
      }
    }
  }
  
  for (auto&& refTkMatches : refTkMap)
  {
    if (refTkMatches.second.size() < 2) // no duplicates
    {
      auto& extra(evt_extras[refTkMatches.second[0]]);
      extra.setDuplicateInfo(0,bool(false));
    }
    else // sort duplicates (ghosts) to keep best one --> most hits, lowest chi2
    {  
      // really should sort on indices with a reduced data structure... this is a hack way to do this for now...
      // e.g. std::tuple<int, int, float>, (label, nHits, chi2)
      TrackVec tmpMatches;
      for (auto&& label : refTkMatches.second) // loop over vector of reco track labels, push back the track with each label 
      {
	tmpMatches.emplace_back(evt_tracks[label]);
      }
      std::sort(tmpMatches.begin(), tmpMatches.end(), sortByHitsChi2); // sort the tracks
      for (auto itrack = 0; itrack < tmpMatches.size(); itrack++) // loop over sorted tracks, now make the vector of sorted labels match
      {
	refTkMatches.second[itrack] = tmpMatches[itrack].label();
      }
      
      int duplicateID = 0;
      for (auto&& label : refTkMatches.second) // loop over vector of reco tracsk 
      {
        auto& extra(evt_extras[label]);
        extra.setDuplicateInfo(duplicateID,bool(true));
        duplicateID++; // used in fake rate trees!
      } 
    }
  }
}

void TTreeValidation::makeSeedTkToRecoTkMaps(Event& ev)
{
  std::lock_guard<std::mutex> locker(glock_); 
  // map seed to reco tracks --> seed track collection assumed to map to itself, unless we make some cuts
  TTreeValidation::mapSeedTkToRecoTk(ev.candidateTracks_,ev.candidateTracksExtra_,seedToBuildMap_);
  TTreeValidation::mapSeedTkToRecoTk(ev.fitTracks_,ev.fitTracksExtra_,seedToFitMap_);
}

void TTreeValidation::mapSeedTkToRecoTk(const TrackVec& evt_tracks, const TrackExtraVec& evt_extras, TkIDToTkIDMap& seedTkMap)
{
  for (auto&& track : evt_tracks)
  {
    seedTkMap[evt_extras[track.label()].seedID()] = track.label();
  }
}

void TTreeValidation::makeRecoTkToRecoTkMaps(Event& ev)
{
  std::lock_guard<std::mutex> locker(glock_);
  TTreeValidation::makeRecoTkToRecoTkMap(buildToFitMap_,ev.candidateTracks_,ev.candidateTracksExtra_,ev.fitTracks_,ev.fitTracksExtra_);
  TTreeValidation::makeRecoTkToRecoTkMap(fitToBuildMap_,ev.fitTracks_,ev.fitTracksExtra_,ev.candidateTracks_,ev.candidateTracksExtra_);
}

void TTreeValidation::makeRecoTkToRecoTkMap(TkIDToTkIDMap& refToPairMap, const TrackVec& reftracks, const TrackExtraVec& refextras, 
					    const TrackVec& pairtracks, const TrackExtraVec& pairextras)
{
  // at this point in the code, the labels of the tracks point their position inside the vector
  // while the seedID is the label prior to relabeling (in reality, this is the MC track ID)
  for (auto&& reftrack : reftracks)
  {
    const auto & refextra = refextras[reftrack.label()];
    for (auto&& pairtrack : pairtracks)
    {
      const auto & pairextra = pairextras[pairtrack.label()];
      if (refextra.seedID() == pairextra.seedID())
      {
	refToPairMap[reftrack.label()] = pairtrack.label(); 
	break;
      }
    }
  }
}

void TTreeValidation::makeCMSSWTkToRecoTksMaps(Event& ev)
{
  std::lock_guard<std::mutex> locker(glock_);
  // can reuse this function
  TTreeValidation::mapRefTkToRecoTks(ev.candidateTracks_,ev.candidateTracksExtra_,cmsswToBuildMap_);
  TTreeValidation::mapRefTkToRecoTks(ev.fitTracks_,ev.fitTracksExtra_,cmsswToFitMap_);
}

void TTreeValidation::makeSeedTkToCMSSWTkMap(Event& ev)
{
  const auto& seedtracks  = ev.seedTracks_;
  const auto& cmsswtracks = ev.cmsswTracks_;
  for (int itrack = 0; itrack < seedtracks.size(); itrack++)
  {
    for (auto&& cmsswtrack : cmsswtracks)
    {
      if (cmsswtrack.label() == itrack) seedToCmsswMap_[seedtracks[itrack].label()] = cmsswtrack.label();
    }
  }
}

void TTreeValidation::storeSeedAndMCID(Event& ev)
{
  const auto& buildtracks = ev.candidateTracks_;
        auto& buildextras = ev.candidateTracksExtra_;

  const auto& fittracks = ev.fitTracks_;
        auto& fitextras = ev.fitTracksExtra_;

  const auto& cmsswtracks = ev.cmsswTracks_;
        auto& cmsswextras = ev.cmsswTracksExtra_;
  
  // first set candidate tracks, use as base for fittracks
  int newlabel = -1;
  for (int itrack = 0; itrack < buildtracks.size(); itrack++)
  {
    auto& extra = buildextras[itrack];
    const int seedID = extra.seedID();

    extra.setmcTrackID(seedID);

    if (seedToCmsswMap_.count(seedID))
    {
      extra.setseedID(seedToCmsswMap_[seedID]);
      if (Config::cmsswMatching == labelBased)
      {
	for (int ctrack = 0; ctrack < cmsswextras.size(); ctrack++)
        {
	  if (cmsswextras[ctrack].seedID() == extra.seedID())
	  {
	    buildToCmsswMap_[itrack] = cmsswtracks[ctrack].label(); // cmsstracks[ctrack].label() == ctrack!
	    break;
	  }
	}
      }
    }
    else 
    {
      extra.setseedID(--newlabel);
    }
  }

  // set according to candidate tracks for fit tracks through map
  for (int itrack = 0; itrack < fittracks.size(); itrack++)
  {
    auto& extra = fitextras[itrack];

    extra.setmcTrackID(buildextras[fitToBuildMap_[itrack]].mcTrackID());
    extra.setseedID   (buildextras[fitToBuildMap_[itrack]].seedID());
  }
}

void TTreeValidation::setupCMSSWMatching(const Event & ev, RedTrackVec & reducedCMSSW, LayIdxIDVecMapMap & cmsswHitIDMap)
{
  // get the cmssw tracks
  const auto& cmsswtracks = ev.cmsswTracks_;
        auto& cmsswextras = ev.cmsswTracksExtra_;

  // resize accordingly
  reducedCMSSW.resize(cmsswtracks.size());

  for (int itrack = 0; itrack < cmsswtracks.size(); itrack++)
  {
    const auto & cmsswtrack = cmsswtracks[itrack];
    const int seedID = cmsswextras[itrack].seedID();
    const SVector6 & params = cmsswtrack.parameters();
    SVector2 tmpv(params[3],params[5]);
    
    HitLayerMap tmpmap;
    for (int ihit = 0; ihit < cmsswtrack.nTotalHits(); ihit++)
    {
      const int lyr = cmsswtrack.getHitLyr(ihit);
      const int idx = cmsswtrack.getHitIdx(ihit);
      
      if (lyr >= 0 && idx >= 0) 
      {
	tmpmap[lyr].push_back(idx);
	cmsswHitIDMap[lyr][idx].push_back(cmsswtrack.label());
      }
    }

    // index inside object is label (as cmsswtracks are now aligned)
    reducedCMSSW[itrack] = ReducedTrack(cmsswtrack.label(),seedID,tmpv,cmsswtrack.momPhi(),tmpmap);
  }
}

int TTreeValidation::getLastFoundHit(const int trackMCHitID, const int mcTrackID, const Event& ev)
{
  int mcHitID = -1;
  if (ev.simHitsInfo_[trackMCHitID].mcTrackID() == mcTrackID)
  {
    mcHitID = trackMCHitID;
  }
  else
  {
    mcHitID = ev.simTracks_[mcTrackID].getMCHitIDFromLayer(ev.layerHits_,ev.simHitsInfo_[trackMCHitID].layer());
  }
  return mcHitID;
}

void TTreeValidation::resetFitBranches()
{
  for(int ilayer = 0; ilayer < Config::nTotalLayers; ++ilayer)
  {
    z_prop_fit_[ilayer]  = -1000.f;
    ez_prop_fit_[ilayer] = -1000.f;
    z_hit_fit_[ilayer]   = -1000.f;
    ez_hit_fit_[ilayer]  = -1000.f;
    z_sim_fit_[ilayer]   = -1000.f;
    ez_sim_fit_[ilayer]  = -1000.f;
    
    pphi_prop_fit_[ilayer]  = -1000.f;
    epphi_prop_fit_[ilayer] = -1000.f;
    pphi_hit_fit_[ilayer]   = -1000.f;
    epphi_hit_fit_[ilayer]  = -1000.f;
    pphi_sim_fit_[ilayer]   = -1000.f;
    epphi_sim_fit_[ilayer]  = -1000.f;

    pt_up_fit_[ilayer]   = -1000.f;
    ept_up_fit_[ilayer]  = -1000.f;
    pt_sim_fit_[ilayer]  = -1000.f;
    ept_sim_fit_[ilayer] = -1000.f;

    mphi_up_fit_[ilayer]   = -1000.f;
    emphi_up_fit_[ilayer]  = -1000.f;
    mphi_sim_fit_[ilayer]  = -1000.f;
    emphi_sim_fit_[ilayer] = -1000.f;

    meta_up_fit_[ilayer]   = -1000.f;
    emeta_up_fit_[ilayer]  = -1000.f;
    meta_sim_fit_[ilayer]  = -1000.f;
    emeta_sim_fit_[ilayer] = -1000.f;
  }  
}

void TTreeValidation::fillFitTree(const Event& ev)
{
  std::lock_guard<std::mutex> locker(glock_); 

  evtid_fit_ = ev.evtID();
  const auto& simtracks = ev.simTracks_;
  const auto& layerhits = ev.layerHits_;
  const auto& simtrackstates = ev.simTrackStates_;
  
  for(auto&& fitvalmapmap : fitValTkMapMap_)
  {
    TTreeValidation::resetFitBranches();
    
    tkid_fit_ = fitvalmapmap.first; // seed id (label) is the same as the mcID
    
    const auto& simtrack = simtracks[tkid_fit_];
    const auto& fitvalmap = fitvalmapmap.second;
    for(int ilayer = 0; ilayer < Config::nTotalLayers; ++ilayer)
    {
      if (fitvalmap.count(ilayer))
      {
	const auto& hit    = layerhits[ilayer][simtrack.getHitIdx(ilayer)];
	const auto& initTS = simtrackstates.at(hit.mcHitID());
	const auto& fitval = fitvalmap.at(ilayer);
	
	z_hit_fit_[ilayer]   = hit.z();
	ez_hit_fit_[ilayer]  = std::sqrt(hit.ezz());
	z_sim_fit_[ilayer]   = initTS.z();
	ez_sim_fit_[ilayer]  = initTS.ezz();
	z_prop_fit_[ilayer]  = fitval.ppz;
	ez_prop_fit_[ilayer] = fitval.eppz;

	pphi_hit_fit_[ilayer]   = hit.phi();
	epphi_hit_fit_[ilayer]  = std::sqrt(hit.ephi());
	pphi_sim_fit_[ilayer]   = initTS.posPhi();
	epphi_sim_fit_[ilayer]  = initTS.eposPhi();
	pphi_prop_fit_[ilayer]  = fitval.ppphi;
	epphi_prop_fit_[ilayer] = fitval.eppphi;
	
	pt_up_fit_[ilayer]   = fitval.upt;
	ept_up_fit_[ilayer]  = fitval.eupt;
	pt_sim_fit_[ilayer]  = initTS.pT();
	ept_sim_fit_[ilayer] = initTS.epT();

	mphi_up_fit_[ilayer]   = fitval.umphi;
	emphi_up_fit_[ilayer]  = fitval.eumphi;
	mphi_sim_fit_[ilayer]  = initTS.momPhi();
	emphi_sim_fit_[ilayer] = initTS.emomPhi();

	meta_up_fit_[ilayer]   = fitval.umeta;
	emeta_up_fit_[ilayer]  = fitval.eumeta;
	meta_sim_fit_[ilayer]  = initTS.momEta();
	emeta_sim_fit_[ilayer] = initTS.emomEta();	
      }
    }
    fittree_->Fill();
  }
}

void TTreeValidation::fillEfficiencyTree(const Event& ev)
{
  std::lock_guard<std::mutex> locker(glock_);

  auto ievt = ev.evtID();
  auto& evt_sim_tracks   = ev.simTracks_;
  auto& evt_seed_tracks  = ev.seedTracks_;
  auto& evt_seed_extras  = ev.seedTracksExtra_;
  auto& evt_build_tracks = ev.candidateTracks_;
  auto& evt_build_extras = ev.candidateTracksExtra_;
  auto& evt_fit_tracks   = ev.fitTracks_;
  auto& evt_fit_extras   = ev.fitTracksExtra_;
  auto& evt_layer_hits   = ev.layerHits_;
  const auto& evt_sim_trackstates = ev.simTrackStates_;

  for (auto&& simtrack : evt_sim_tracks)
  {
    evtID_eff_ = ievt;
    mcID_eff_  = simtrack.label();

    // generated values
    x_mc_gen_eff_ = simtrack.x();
    y_mc_gen_eff_ = simtrack.y();
    z_mc_gen_eff_ = simtrack.z();

    pt_mc_gen_eff_  = simtrack.pT(); 
    phi_mc_gen_eff_ = simtrack.momPhi();
    eta_mc_gen_eff_ = simtrack.momEta();
    nHits_mc_eff_   = simtrack.nFoundHits(); // could be that the sim track skips layers!
    lastlyr_mc_eff_ = simtrack.getLastFoundHitLyr();

    // matched seed track
    if (simToSeedMap_.count(mcID_eff_) && simtrack.isFindable()) // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToSeedMap_[matched SimID][first element in vector]
    {
      auto& seedtrack = evt_seed_tracks[simToSeedMap_[mcID_eff_][0]]; // returns seedTrack best matched to sim track
      auto& seedextra = evt_seed_extras[seedtrack.label()]; // returns track extra best aligned with seed track
      mcmask_seed_eff_ = 1; // quick logic for matched

      seedID_seed_eff_ = seedextra.seedID(); 

      // use this to access correct sim track layer params
      const int mcHitID = TTreeValidation::getLastFoundHit(seedtrack.getLastFoundMCHitID(evt_layer_hits),mcID_eff_,ev);
      if (mcHitID >= 0 && Config::readSimTrackStates)
      {
	const TrackState & initLayTS = evt_sim_trackstates[mcHitID];

	pt_mc_seed_eff_  = initLayTS.pT();
	phi_mc_seed_eff_ = initLayTS.momPhi();
	eta_mc_seed_eff_ = initLayTS.momEta();
	helixchi2_seed_eff_ = computeHelixChi2(initLayTS.parameters,seedtrack.parameters(),seedtrack.errors());
      }	
      else
      {
	pt_mc_seed_eff_  = -101;
	phi_mc_seed_eff_ = -101;
	eta_mc_seed_eff_ = -101;
	helixchi2_seed_eff_ = -101;
      }

      // last hit info
      const Hit& lasthit = evt_layer_hits[seedtrack.getLastFoundHitLyr()][seedtrack.getLastFoundHitIdx()];
      xhit_seed_eff_ = lasthit.x(); 
      yhit_seed_eff_ = lasthit.y(); 
      zhit_seed_eff_ = lasthit.z(); 

      pt_seed_eff_   = seedtrack.pT();
      ept_seed_eff_  = seedtrack.epT();
      phi_seed_eff_  = seedtrack.momPhi();
      ephi_seed_eff_ = seedtrack.emomPhi();
      eta_seed_eff_  = seedtrack.momEta();
      eeta_seed_eff_ = seedtrack.emomEta();

      // rest of mc info
      nHits_seed_eff_           = seedtrack.nFoundHits();
      nHitsMatched_seed_eff_    = seedextra.nHitsMatched();
      fracHitsMatched_seed_eff_ = seedextra.fracHitsMatched();
      lastlyr_seed_eff_         = seedtrack.getLastFoundHitLyr();

      hitchi2_seed_eff_ = seedtrack.chi2(); // currently not being used

      duplmask_seed_eff_   = seedextra.isDuplicate(); 
      nTkMatches_seed_eff_ = simToSeedMap_[mcID_eff_].size(); // n reco matches to this sim track.
    }
    else // unmatched simTracks ... put -99 for all reco values to denote unmatched
    {
      mcmask_seed_eff_ = (simtrack.isFindable() ? 0 : -1); // quick logic for not matched
      
      seedID_seed_eff_ = -99;
      
      pt_mc_seed_eff_  = -99;
      phi_mc_seed_eff_ = -99;
      eta_mc_seed_eff_ = -99;
      helixchi2_seed_eff_ = -99;

      xhit_seed_eff_ = -2000;
      yhit_seed_eff_ = -2000;
      zhit_seed_eff_ = -2000;

      pt_seed_eff_   = -99;
      ept_seed_eff_  = -99;
      phi_seed_eff_  = -99;
      ephi_seed_eff_ = -99;
      eta_seed_eff_  = -99;
      eeta_seed_eff_ = -99;

      nHits_seed_eff_           = -99;
      nHitsMatched_seed_eff_    = -99;
      fracHitsMatched_seed_eff_ = -99;
      lastlyr_seed_eff_         = -99;
 
      hitchi2_seed_eff_   = -99;
      
      duplmask_seed_eff_   = -1; // mask means unmatched sim track
      nTkMatches_seed_eff_ = -99; // unmatched
    }

    // matched build track
    if (simToBuildMap_.count(mcID_eff_) && simtrack.isFindable()) // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToBuildMap_[matched SimID][first element in vector]
    {
      auto& buildtrack = evt_build_tracks[simToBuildMap_[mcID_eff_][0]]; // returns buildTrack best matched to sim track
      auto& buildextra = evt_build_extras[buildtrack.label()]; // returns track extra best aligned with build track
      mcmask_build_eff_ = 1; // quick logic for matched

      seedID_build_eff_ = buildextra.seedID(); 

      // use this to access correct sim track layer params
      const int mcHitID = TTreeValidation::getLastFoundHit(buildtrack.getLastFoundMCHitID(evt_layer_hits),mcID_eff_,ev);
      if (mcHitID >= 0 && Config::readSimTrackStates)
      {
	const TrackState & initLayTS = evt_sim_trackstates[mcHitID];

	pt_mc_build_eff_  = initLayTS.pT();
	phi_mc_build_eff_ = initLayTS.momPhi();
	eta_mc_build_eff_ = initLayTS.momEta();
	helixchi2_build_eff_ = computeHelixChi2(initLayTS.parameters,buildtrack.parameters(),buildtrack.errors());
      }	
      else
      {
	pt_mc_build_eff_  = -101;
	phi_mc_build_eff_ = -101;
	eta_mc_build_eff_ = -101;
	helixchi2_build_eff_ = -101;
      }

      // last hit info
      const Hit& lasthit = evt_layer_hits[buildtrack.getLastFoundHitLyr()][buildtrack.getLastFoundHitIdx()];
      xhit_build_eff_ = lasthit.x(); 
      yhit_build_eff_ = lasthit.y(); 
      zhit_build_eff_ = lasthit.z(); 

      pt_build_eff_   = buildtrack.pT();
      ept_build_eff_  = buildtrack.epT();
      phi_build_eff_  = buildtrack.momPhi();
      ephi_build_eff_ = buildtrack.emomPhi();
      eta_build_eff_  = buildtrack.momEta();
      eeta_build_eff_ = buildtrack.emomEta();
      
      nHits_build_eff_           = buildtrack.nFoundHits();
      nHitsMatched_build_eff_    = buildextra.nHitsMatched();
      fracHitsMatched_build_eff_ = buildextra.fracHitsMatched();
      lastlyr_build_eff_         = buildtrack.getLastFoundHitLyr();

      hitchi2_build_eff_ = buildtrack.chi2(); 

      duplmask_build_eff_   = buildextra.isDuplicate(); 
      nTkMatches_build_eff_ = simToBuildMap_[mcID_eff_].size(); // n reco matches to this sim track.
    }
    else // unmatched simTracks ... put -99 for all reco values to denote unmatched
    {
      mcmask_build_eff_ = (simtrack.isFindable() ? 0 : -1); // quick logic for not matched

      seedID_build_eff_ = -99;

      pt_mc_build_eff_  = -99;
      phi_mc_build_eff_ = -99;
      eta_mc_build_eff_ = -99;
      helixchi2_build_eff_ = -99;

      xhit_build_eff_ = -2000;
      yhit_build_eff_ = -2000;
      zhit_build_eff_ = -2000;

      pt_build_eff_   = -99;
      ept_build_eff_  = -99;
      phi_build_eff_  = -99;
      ephi_build_eff_ = -99;
      eta_build_eff_  = -99;
      eeta_build_eff_ = -99;

      nHits_build_eff_           = -99;
      nHitsMatched_build_eff_    = -99;
      fracHitsMatched_build_eff_ = -99;
      lastlyr_build_eff_         = -99;

      hitchi2_build_eff_   = -99;
      
      duplmask_build_eff_   = -1; // mask means unmatched sim track
      nTkMatches_build_eff_ = -99; // unmatched
    }
    
    // matched fit track
    if (simToFitMap_.count(mcID_eff_) && simtrack.isFindable()) // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToFitMap_[matched SimID][first element in vector]
    {
      auto& fittrack = evt_fit_tracks[simToFitMap_[mcID_eff_][0]]; // returns fitTrack best matched to sim track
      auto& fitextra = evt_fit_extras[fittrack.label()]; // returns track extra best aligned with fit track
      mcmask_fit_eff_ = 1; // quick logic for matched

      seedID_fit_eff_ = fitextra.seedID(); 

      // use this to access correct sim track layer params
      const int mcHitID = TTreeValidation::getLastFoundHit(fittrack.getLastFoundMCHitID(evt_layer_hits),mcID_eff_,ev);
      if (mcHitID >= 0 && Config::readSimTrackStates)
      {
	const TrackState & initLayTS = evt_sim_trackstates[mcHitID];

	pt_mc_fit_eff_  = initLayTS.pT();
	phi_mc_fit_eff_ = initLayTS.momPhi();
	eta_mc_fit_eff_ = initLayTS.momEta();
	helixchi2_fit_eff_ = computeHelixChi2(initLayTS.parameters,fittrack.parameters(),fittrack.errors());
      }	
      else
      {
	pt_mc_fit_eff_  = -101;
	phi_mc_fit_eff_ = -101;
	eta_mc_fit_eff_ = -101;
	helixchi2_fit_eff_ = -101;
      }
      
      // last hit info
      const Hit& lasthit = evt_layer_hits[fittrack.getLastFoundHitLyr()][fittrack.getLastFoundHitIdx()];
      xhit_fit_eff_ = lasthit.x(); 
      yhit_fit_eff_ = lasthit.y(); 
      zhit_fit_eff_ = lasthit.z(); 

      pt_fit_eff_   = fittrack.pT();
      ept_fit_eff_  = fittrack.epT();
      phi_fit_eff_  = fittrack.momPhi();
      ephi_fit_eff_ = fittrack.emomPhi();
      eta_fit_eff_  = fittrack.momEta();
      eeta_fit_eff_ = fittrack.emomEta();
      
      // rest of mc info
      nHits_fit_eff_           = fittrack.nFoundHits();
      nHitsMatched_fit_eff_    = fitextra.nHitsMatched();
      fracHitsMatched_fit_eff_ = fitextra.fracHitsMatched();
      lastlyr_fit_eff_         = fittrack.getLastFoundHitLyr();

      hitchi2_fit_eff_ = -10; //fittrack.chi2(); // currently not being used

      duplmask_fit_eff_   = fitextra.isDuplicate(); 
      nTkMatches_fit_eff_ = simToFitMap_[mcID_eff_].size(); // n reco matches to this sim track.
    }
    else // unmatched simTracks ... put -99 for all reco values to denote unmatched
    {
      mcmask_fit_eff_ = (simtrack.isFindable() ? 0 : -1); // quick logic for not matched

      seedID_fit_eff_ = -99;

      pt_mc_fit_eff_  = -99;
      phi_mc_fit_eff_ = -99;
      eta_mc_fit_eff_ = -99;
      helixchi2_fit_eff_ = -99;

      xhit_fit_eff_ = -2000;
      yhit_fit_eff_ = -2000;
      zhit_fit_eff_ = -2000;

      pt_fit_eff_   = -99;
      ept_fit_eff_  = -99;
      phi_fit_eff_  = -99;
      ephi_fit_eff_ = -99;
      eta_fit_eff_  = -99;
      eeta_fit_eff_ = -99;

      nHits_fit_eff_           = -99;
      nHitsMatched_fit_eff_    = -99;
      fracHitsMatched_fit_eff_ = -99;
      lastlyr_fit_eff_         = -99;

      hitchi2_fit_eff_   = -99;
      
      duplmask_fit_eff_   = -1; // mask means unmatched sim track
      nTkMatches_fit_eff_ = -99; // unmatched
    }

    efftree_->Fill(); // fill it once per sim track!
  }
}

void TTreeValidation::fillFakeRateTree(const Event& ev)
{
  std::lock_guard<std::mutex> locker(glock_);

  auto ievt = ev.evtID();
  auto& evt_sim_tracks   = ev.simTracks_; // store sim info at that final layer!!! --> gen info stored only in eff tree
  auto& evt_seed_tracks  = ev.seedTracks_;
  auto& evt_seed_extras  = ev.seedTracksExtra_;
  auto& evt_build_tracks = ev.candidateTracks_;
  auto& evt_build_extras = ev.candidateTracksExtra_;
  auto& evt_fit_tracks   = ev.fitTracks_;
  auto& evt_fit_extras   = ev.fitTracksExtra_;
  auto& evt_layer_hits   = ev.layerHits_;
  auto& evt_sim_trackstates = ev.simTrackStates_;

  for (auto&& seedtrack : evt_seed_tracks)
  { 
    evtID_FR_       = ievt;
    auto& seedextra = evt_seed_extras[seedtrack.label()];
    seedID_FR_      = seedextra.seedID();

    // seed info
    seedmask_seed_FR_ = 1; // automatically set to 1, because at the moment no cuts on seeds after conformal+KF fit.  seed triplets filtered by RZ chi2 before fitting. 

    // last hit info
    const Hit& lasthit = evt_layer_hits[seedtrack.getLastFoundHitLyr()][seedtrack.getLastFoundHitIdx()];
    xhit_seed_FR_ = lasthit.x(); 
    yhit_seed_FR_ = lasthit.y(); 
    zhit_seed_FR_ = lasthit.z(); 
    
    pt_seed_FR_   = seedtrack.pT();
    ept_seed_FR_  = seedtrack.epT();
    phi_seed_FR_  = seedtrack.momPhi();
    ephi_seed_FR_ = seedtrack.emomPhi();
    eta_seed_FR_  = seedtrack.momEta();
    eeta_seed_FR_ = seedtrack.emomEta();

    nHits_seed_FR_           = seedtrack.nFoundHits();
    nHitsMatched_seed_FR_    = seedextra.nHitsMatched();
    fracHitsMatched_seed_FR_ = seedextra.fracHitsMatched();
    lastlyr_seed_FR_         = seedtrack.getLastFoundHitLyr();

    hitchi2_seed_FR_ = seedtrack.chi2(); //--> not currently used

    // sim info for seed track
    mcID_seed_FR_ = seedextra.mcTrackID();
    if (mcID_seed_FR_ >= 0) // seed track matched to seed and sim 
    {
      mcmask_seed_FR_ = 1; // matched track to sim
    }
    else 
    {
      if (Config::inclusiveShorts) 
      {
	if      (mcID_seed_FR_ ==  -1 || mcID_seed_FR_ ==  -5 || 
		 mcID_seed_FR_ ==  -8 || mcID_seed_FR_ ==  -9)  
	{
	  mcmask_seed_FR_ = 0;
	}
	else if (mcID_seed_FR_ ==  -2 || mcID_seed_FR_ == -10 ||
		 mcID_seed_FR_ == -11)
	{
	  mcmask_seed_FR_ = 2; 
	}
	else // mcID == -3,-4,-6,-7,-12,-13
	{
	  mcmask_seed_FR_ = -1;
	}
      }
      else // only count long tracks
      {
	if      (mcID_seed_FR_ == -1 || mcID_seed_FR_ == -9) 
	{
	  mcmask_seed_FR_ = 0;
	}
	else // mcID == -2,-3,-5,-6,-7,-8,-10,-11,-12,-13
	{
	  mcmask_seed_FR_ = -1; 
	}
      }
    } // end check over not matched
    
    if (mcmask_seed_FR_ == 1) // matched track to sim
    {
      auto& simtrack = evt_sim_tracks[mcID_seed_FR_];

      const int mcHitID = TTreeValidation::getLastFoundHit(seedtrack.getLastFoundMCHitID(evt_layer_hits),mcID_seed_FR_,ev);
      if (mcHitID >= 0 && Config::readSimTrackStates)
      {
	const TrackState & initLayTS = evt_sim_trackstates[mcHitID];
	pt_mc_seed_FR_  = initLayTS.pT();
	phi_mc_seed_FR_ = initLayTS.momPhi();
	eta_mc_seed_FR_ = initLayTS.momEta();
	helixchi2_seed_FR_ = computeHelixChi2(initLayTS.parameters,seedtrack.parameters(),seedtrack.errors());
      }	
      else
      {
	pt_mc_seed_FR_  = -101;
	phi_mc_seed_FR_ = -101;
	eta_mc_seed_FR_ = -101;
	helixchi2_seed_FR_ = -101;
      }

      nHits_mc_seed_FR_   = simtrack.nFoundHits();

      lastlyr_mc_seed_FR_ = simtrack.getLastFoundHitLyr();

      duplmask_seed_FR_   = seedextra.isDuplicate();
      iTkMatches_seed_FR_ = seedextra.duplicateID(); // ith duplicate seed track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"      
    }
    else
    {
      // -99 for all sim info for reco tracks not associated to reco tracks
      pt_mc_seed_FR_  = -99;
      phi_mc_seed_FR_ = -99;
      eta_mc_seed_FR_ = -99;
      helixchi2_seed_FR_ = -99;
      
      nHits_mc_seed_FR_ = -99;
      lastlyr_mc_seed_FR_ = -99;

      duplmask_seed_FR_   = -1;
      iTkMatches_seed_FR_ = -99;  
    }

    //==========================//

    // fill build information if track still alive
    if (seedToBuildMap_.count(seedID_FR_))
    {
      seedmask_build_FR_ = 1; // quick logic

      auto& buildtrack = evt_build_tracks[seedToBuildMap_[seedID_FR_]];
      auto& buildextra = evt_build_extras[buildtrack.label()];

      // last hit info
      const Hit& lasthit = evt_layer_hits[buildtrack.getLastFoundHitLyr()][buildtrack.getLastFoundHitIdx()];
      xhit_build_FR_ = lasthit.x(); 
      yhit_build_FR_ = lasthit.y(); 
      zhit_build_FR_ = lasthit.z(); 

      pt_build_FR_   = buildtrack.pT();
      ept_build_FR_  = buildtrack.epT();
      phi_build_FR_  = buildtrack.momPhi();
      ephi_build_FR_ = buildtrack.emomPhi();
      eta_build_FR_  = buildtrack.momEta();
      eeta_build_FR_ = buildtrack.emomEta();

      nHits_build_FR_           = buildtrack.nFoundHits();
      nHitsMatched_build_FR_    = buildextra.nHitsMatched();
      fracHitsMatched_build_FR_ = buildextra.fracHitsMatched();
      lastlyr_build_FR_         = buildtrack.getLastFoundHitLyr();

      hitchi2_build_FR_ = buildtrack.chi2();

      // sim info for build track
      mcID_build_FR_  = buildextra.mcTrackID();
      if (mcID_build_FR_ >= 0) // build track matched to seed and sim 
      {
	mcmask_build_FR_ = 1; // matched track to sim
      }
      else 
      {
	if (Config::inclusiveShorts) 
        {
	  if      (mcID_build_FR_ ==  -1 || mcID_build_FR_ ==  -5 ||
		   mcID_build_FR_ ==  -8 || mcID_build_FR_ ==  -9)  
	  {
	    mcmask_build_FR_ = 0;
	  }
	  else if (mcID_build_FR_ ==  -2 || mcID_build_FR_ == -10 ||
		   mcID_build_FR_ == -11)
	  {
	    mcmask_build_FR_ = 2; 
	  }
	  else // mcID == -3,-4,-6,-7,-12,-13
	  {
	    mcmask_build_FR_ = -1;
	  }
	}
	else // only count long tracks
        {
	  if      (mcID_build_FR_ == -1 || mcID_build_FR_ == -9) 
	  {
	    mcmask_build_FR_ = 0;
	  }
	  else // mcID == -2,-3,-4,-5,-6,-7,-8,-10,-11,-12,-13
	  {
	    mcmask_build_FR_ = -1; 
	  }
	}
      } // end check over not matched

      if (mcmask_build_FR_ == 1) // build track matched to seed and sim 
      {
	auto& simtrack = evt_sim_tracks[mcID_build_FR_];

	const int mcHitID = TTreeValidation::getLastFoundHit(buildtrack.getLastFoundMCHitID(evt_layer_hits),mcID_build_FR_,ev);
	if (mcHitID >= 0 && Config::readSimTrackStates)
        {
	  const TrackState & initLayTS = evt_sim_trackstates[mcHitID];
	  pt_mc_build_FR_  = initLayTS.pT();
	  phi_mc_build_FR_ = initLayTS.momPhi();
	  eta_mc_build_FR_ = initLayTS.momEta();
	  helixchi2_build_FR_ = computeHelixChi2(initLayTS.parameters,buildtrack.parameters(),buildtrack.errors());
	}	
	else 
        {
	  pt_mc_build_FR_  = -101;
	  phi_mc_build_FR_ = -101;
	  eta_mc_build_FR_ = -101;
	  helixchi2_build_FR_ = -101;
	}

	nHits_mc_build_FR_   = simtrack.nFoundHits();
	lastlyr_mc_build_FR_ = simtrack.getLastFoundHitLyr();

	duplmask_build_FR_   = buildextra.isDuplicate();
	iTkMatches_build_FR_ = buildextra.duplicateID(); // ith duplicate build track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"      
      }
      else // build track matched only to seed not to sim
      {
	// -99 for all sim info for reco tracks not associated to reco tracks
	pt_mc_build_FR_  = -99;
	phi_mc_build_FR_ = -99;
	eta_mc_build_FR_ = -99;
	helixchi2_build_FR_ = -99;

	nHits_mc_build_FR_ = -99;
	lastlyr_mc_build_FR_ = -99;

	duplmask_build_FR_   = -1;
	iTkMatches_build_FR_ = -99;
      } // matched seed to build, not build to sim
    }

    else  // seed has no matching build track (therefore no matching sim to build track)
    { 
      seedmask_build_FR_ = 0; // quick logic

      // -3000 for position info if no build track for seed
      xhit_build_FR_ = -3000;
      yhit_build_FR_ = -3000;
      zhit_build_FR_ = -3000;

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
      lastlyr_build_FR_ = -100;

      hitchi2_build_FR_  = -100; 
      
      // keep -100 for all sim variables as no such reco exists for this seed
      mcmask_build_FR_ = -2; // do not want to count towards build FR
      mcID_build_FR_   = -100;
	
      pt_mc_build_FR_  = -100;
      phi_mc_build_FR_ = -100;
      eta_mc_build_FR_ = -100;
      helixchi2_build_FR_ = -100;

      nHits_mc_build_FR_ = -100;
      lastlyr_mc_build_FR_ = -100;

      duplmask_build_FR_   = -2;
      iTkMatches_build_FR_ = -100;
    }

    //============================// fit tracks
    if (seedToFitMap_.count(seedID_FR_))
    {
      seedmask_fit_FR_ = 1; // quick logic

      auto& fittrack = evt_fit_tracks[seedToFitMap_[seedID_FR_]];
      auto& fitextra = evt_fit_extras[fittrack.label()];

      // last hit info
      const Hit& lasthit = evt_layer_hits[fittrack.getLastFoundHitLyr()][fittrack.getLastFoundHitIdx()];
      xhit_fit_FR_ = lasthit.x(); 
      yhit_fit_FR_ = lasthit.y(); 
      zhit_fit_FR_ = lasthit.z(); 

      pt_fit_FR_   = fittrack.pT();
      ept_fit_FR_  = fittrack.epT();
      phi_fit_FR_  = fittrack.momPhi();
      ephi_fit_FR_ = fittrack.emomPhi();
      eta_fit_FR_  = fittrack.momEta();
      eeta_fit_FR_ = fittrack.emomEta();

      nHits_fit_FR_           = fittrack.nFoundHits();
      nHitsMatched_fit_FR_    = fitextra.nHitsMatched();
      fracHitsMatched_fit_FR_ = fitextra.fracHitsMatched();
      lastlyr_fit_FR_         = fittrack.getLastFoundHitLyr();

      hitchi2_fit_FR_ = -10; //fittrack.chi2() --> currently not used

      // sim info for fit track
      mcID_fit_FR_  = fitextra.mcTrackID();
      if (mcID_fit_FR_ >= 0) // fit track matched to seed and sim 
      {
	mcmask_fit_FR_ = 1; // matched track to sim
      }
      else 
      {
	if (Config::inclusiveShorts) 
        {
	  if      (mcID_fit_FR_ ==  -1 || mcID_fit_FR_ ==  -5 ||
		   mcID_fit_FR_ ==  -8 || mcID_fit_FR_ ==  -9)  
	  {
	    mcmask_fit_FR_ = 0;
	  }
	  else if (mcID_fit_FR_ ==  -2 || mcID_fit_FR_ == -10 ||
		   mcID_fit_FR_ == -11)
	  {
	    mcmask_fit_FR_ = 2; 
	  }
	  else // mcID == -3,-4,-6,-7,-12,-13
	  {
	    mcmask_fit_FR_ = -1;
	  }
	}
	else // only count long tracks
        {
	  if      (mcID_fit_FR_ == -1 || mcID_fit_FR_ == -9) 
	  {
	    mcmask_fit_FR_ = 0;
	  }
	  else // mcID == -2,-3,-4,-5,-6,-7,-8,-10,-11,-12,-13
	  {
	    mcmask_fit_FR_ = -1; 
	  }
	}
      } // end check over not matched
   
      if (mcmask_fit_FR_ == 1) // fit track matched to seed and sim 
      {
	auto& simtrack = evt_sim_tracks[mcID_fit_FR_];
      
	const int mcHitID = TTreeValidation::getLastFoundHit(fittrack.getLastFoundMCHitID(evt_layer_hits),mcID_fit_FR_,ev); // only works for outward fit for now
	if (mcHitID >= 0 && Config::readSimTrackStates)
        {
	  const TrackState & initLayTS = evt_sim_trackstates[mcHitID];
	  pt_mc_fit_FR_  = initLayTS.pT();
	  phi_mc_fit_FR_ = initLayTS.momPhi();
	  eta_mc_fit_FR_ = initLayTS.momEta();
	  helixchi2_fit_FR_ = computeHelixChi2(initLayTS.parameters,fittrack.parameters(),fittrack.errors());
	}	
	else
        {
	  pt_mc_fit_FR_  = -101;
	  phi_mc_fit_FR_ = -101;
	  eta_mc_fit_FR_ = -101;
	  helixchi2_fit_FR_ = -101;
	}

	nHits_mc_fit_FR_   = simtrack.nFoundHits();
	lastlyr_mc_fit_FR_ = simtrack.getLastFoundHitLyr();

	duplmask_fit_FR_   = fitextra.isDuplicate();
	iTkMatches_fit_FR_ = fitextra.duplicateID(); // ith duplicate fit track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"
      }
      else // fit track matched only to seed not to sim
      {
	// -99 for all sim info for reco tracks not associated to reco tracks
	pt_mc_fit_FR_  = -99;
	phi_mc_fit_FR_ = -99;
	eta_mc_fit_FR_ = -99;
	helixchi2_fit_FR_ = -99;
	
	nHits_mc_fit_FR_ = -99;
	lastlyr_mc_fit_FR_ = -99;

	duplmask_fit_FR_   = -1;
	iTkMatches_fit_FR_ = -99;
      } // matched seed to fit, not fit to sim
    }

    else // seed has no matching fit track (therefore no matching sim to fit track)
    {
      seedmask_fit_FR_ = 0; // quick logic

      // -3000 for position info if no fit track for seed
      xhit_fit_FR_ = -3000;
      yhit_fit_FR_ = -3000;
      zhit_fit_FR_ = -3000;

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
      lastlyr_fit_FR_ = -100;

      hitchi2_fit_FR_  = -100; 
      
      // keep -100 for all sim variables as no such reco exists for this seed
      mcmask_fit_FR_ = -2; // do not want to count towards fit FR
      mcID_fit_FR_   = -100;
	
      pt_mc_fit_FR_  = -100;
      phi_mc_fit_FR_ = -100;
      eta_mc_fit_FR_ = -100;
      helixchi2_fit_FR_ = -100;

      nHits_mc_fit_FR_ = -100;
      lastlyr_mc_fit_FR_ = -100;

      duplmask_fit_FR_   = -2;
      iTkMatches_fit_FR_ = -100;
    }

    frtree_->Fill(); // fill once per seed!
  }// end of seed to seed loop
}

void TTreeValidation::fillConfigTree()
{
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

void TTreeValidation::fillCMSSWEfficiencyTree(const Event& ev)
{
  std::lock_guard<std::mutex> locker(glock_);

  auto ievt = ev.evtID();
  auto& evt_cmssw_tracks = ev.cmsswTracks_;
  auto& evt_cmssw_extras = ev.cmsswTracksExtra_;
  auto& evt_build_tracks = ev.candidateTracks_;
  auto& evt_build_extras = ev.candidateTracksExtra_;
  auto& evt_fit_tracks   = ev.fitTracks_;
  auto& evt_fit_extras   = ev.fitTracksExtra_;
  auto& evt_layer_hits   = ev.layerHits_;

  for (auto&& cmsswtrack : evt_cmssw_tracks)
  {
    const auto& cmsswextra = evt_cmssw_extras[cmsswtrack.label()];
    
    evtID_ceff_        = ievt;
    cmsswID_ceff_      = cmsswtrack.label();
    seedID_cmssw_ceff_ = cmsswextra.seedID();

    // PCA parameters
    x_cmssw_ceff_ = cmsswtrack.x();
    y_cmssw_ceff_ = cmsswtrack.y();
    z_cmssw_ceff_ = cmsswtrack.z();

    pt_cmssw_ceff_  = cmsswtrack.pT(); 
    phi_cmssw_ceff_ = cmsswtrack.momPhi();
    eta_cmssw_ceff_ = cmsswtrack.momEta();

    nHits_cmssw_ceff_   = cmsswtrack.nFoundHits(); 
    nLayers_cmssw_ceff_ = cmsswtrack.nUniqueLayers();
    lastlyr_cmssw_ceff_ = cmsswtrack.getLastFoundHitLyr();

    // matched build track
    if (cmsswToBuildMap_.count(cmsswID_ceff_) && cmsswtrack.isFindable()) // recoToCmssw match : save best match --> most hits, lowest chi2, i.e. cmsswToBuildMap_[matched CmsswID][first element in vector]
    {
      const auto& buildtrack = evt_build_tracks[cmsswToBuildMap_[cmsswID_ceff_][0]]; // returns buildTrack best matched to cmssw track
      const auto& buildextra = evt_build_extras[buildtrack.label()]; // returns track extra best aligned with build track
      cmsswmask_build_ceff_ = 1; // quick logic for matched

      seedID_build_ceff_    = buildextra.seedID(); 
      mcTrackID_build_ceff_ = buildextra.mcTrackID(); 

      // track parameters
      pt_build_ceff_   = buildtrack.pT();
      ept_build_ceff_  = buildtrack.epT();
      phi_build_ceff_  = buildtrack.momPhi();
      ephi_build_ceff_ = buildtrack.emomPhi();
      eta_build_ceff_  = buildtrack.momEta();
      eeta_build_ceff_ = buildtrack.emomEta();

      // hit/layer info
      nHits_build_ceff_           = buildtrack.nFoundHits();
      nLayers_build_ceff_         = buildtrack.nUniqueLayers();
      nHitsMatched_build_ceff_    = buildextra.nHitsMatched();
      fracHitsMatched_build_ceff_ = buildextra.fracHitsMatched();
      lastlyr_build_ceff_         = buildtrack.getLastFoundHitLyr();

      // hit info
      const Hit& lasthit = evt_layer_hits[buildtrack.getLastFoundHitLyr()][buildtrack.getLastFoundHitIdx()];
      xhit_build_ceff_ = lasthit.x(); 
      yhit_build_ceff_ = lasthit.y(); 
      zhit_build_ceff_ = lasthit.z(); 

      // chi2 info
      hitchi2_build_ceff_ = buildtrack.chi2(); 
      helixchi2_build_ceff_ = buildextra.helixChi2();

      // swim dphi
      dphi_build_ceff_ = buildextra.dPhi(); 

      // duplicate info
      duplmask_build_ceff_   = buildextra.isDuplicate(); 
      nTkMatches_build_ceff_ = cmsswToBuildMap_[cmsswID_ceff_].size(); // n reco matches to this cmssw track.
    }
    else // unmatched cmsswtracks ... put -99 for all reco values to denote unmatched
    {
      cmsswmask_build_ceff_ = (cmsswtrack.isFindable() ? 0 : -1); // quick logic for not matched
      
      seedID_build_ceff_    = -99;
      mcTrackID_build_ceff_ = -99;

      pt_build_ceff_   = -99;
      ept_build_ceff_  = -99;
      phi_build_ceff_  = -99;
      ephi_build_ceff_ = -99;
      eta_build_ceff_  = -99;
      eeta_build_ceff_ = -99;

      nHits_build_ceff_           = -99;
      nLayers_build_ceff_         = -99;
      nHitsMatched_build_ceff_    = -99;
      fracHitsMatched_build_ceff_ = -99;
      lastlyr_build_ceff_         = -99;
 
      xhit_build_ceff_ = -2000;
      yhit_build_ceff_ = -2000;
      zhit_build_ceff_ = -2000;

      hitchi2_build_ceff_   = -99;
      helixchi2_build_ceff_ = -99;
      
      dphi_build_ceff_ = -99;

      duplmask_build_ceff_   = -1; // mask means unmatched cmssw track
      nTkMatches_build_ceff_ = -99; // unmatched
    }

    // matched fit track
    if (cmsswToFitMap_.count(cmsswID_ceff_) && cmsswtrack.isFindable()) // recoToCmssw match : save best match --> most hits, lowest chi2, i.e. cmsswToFitMap_[matched CmsswID][first element in vector]
    {
      const auto& fittrack = evt_fit_tracks[cmsswToFitMap_[cmsswID_ceff_][0]]; // returns fitTrack best matched to cmssw track
      const auto& fitextra = evt_fit_extras[fittrack.label()]; // returns track extra best aligned with fit track
      cmsswmask_fit_ceff_ = 1; // quick logic for matched

      seedID_fit_ceff_    = fitextra.seedID(); 
      mcTrackID_fit_ceff_ = fitextra.mcTrackID(); 

      // track parameters
      pt_fit_ceff_   = fittrack.pT();
      ept_fit_ceff_  = fittrack.epT();
      phi_fit_ceff_  = fittrack.momPhi();
      ephi_fit_ceff_ = fittrack.emomPhi();
      eta_fit_ceff_  = fittrack.momEta();
      eeta_fit_ceff_ = fittrack.emomEta();

      // hit/layer info
      nHits_fit_ceff_           = fittrack.nFoundHits();
      nLayers_fit_ceff_         = fittrack.nUniqueLayers();
      nHitsMatched_fit_ceff_    = fitextra.nHitsMatched();
      fracHitsMatched_fit_ceff_ = fitextra.fracHitsMatched();
      lastlyr_fit_ceff_         = fittrack.getLastFoundHitLyr();

      // hit info
      const Hit& lasthit = evt_layer_hits[fittrack.getLastFoundHitLyr()][fittrack.getLastFoundHitIdx()];
      xhit_fit_ceff_ = lasthit.x(); 
      yhit_fit_ceff_ = lasthit.y(); 
      zhit_fit_ceff_ = lasthit.z(); 

      // chi2 info
      hitchi2_fit_ceff_ = fittrack.chi2(); 
      helixchi2_fit_ceff_ = fitextra.helixChi2();

      // swim dphi
      dphi_fit_ceff_ = fitextra.dPhi(); 

      // duplicate info
      duplmask_fit_ceff_   = fitextra.isDuplicate(); 
      nTkMatches_fit_ceff_ = cmsswToFitMap_[cmsswID_ceff_].size(); // n reco matches to this cmssw track.
    }
    else // unmatched cmsswtracks ... put -99 for all reco values to denote unmatched
    {
      cmsswmask_fit_ceff_ = (cmsswtrack.isFindable() ? 0 : -1); // quick logic for not matched
      
      seedID_fit_ceff_    = -99;
      mcTrackID_fit_ceff_ = -99;

      pt_fit_ceff_   = -99;
      ept_fit_ceff_  = -99;
      phi_fit_ceff_  = -99;
      ephi_fit_ceff_ = -99;
      eta_fit_ceff_  = -99;
      eeta_fit_ceff_ = -99;

      nHits_fit_ceff_           = -99;
      nLayers_fit_ceff_         = -99;
      nHitsMatched_fit_ceff_    = -99;
      fracHitsMatched_fit_ceff_ = -99;
      lastlyr_fit_ceff_         = -99;
 
      xhit_fit_ceff_ = -2000;
      yhit_fit_ceff_ = -2000;
      zhit_fit_ceff_ = -2000;

      hitchi2_fit_ceff_   = -99;
      helixchi2_fit_ceff_ = -99;
      
      dphi_fit_ceff_ = -99;

      duplmask_fit_ceff_   = -1; // mask means unmatched cmssw track
      nTkMatches_fit_ceff_ = -99; // unmatched
    }

    cmsswefftree_->Fill();
  } 
}

void TTreeValidation::fillCMSSWFakeRateTree(const Event& ev)
{
  std::lock_guard<std::mutex> locker(glock_);

  auto ievt = ev.evtID();
  auto& evt_cmssw_tracks = ev.cmsswTracks_;
  auto& evt_cmssw_extras = ev.cmsswTracksExtra_;
  auto& evt_build_tracks = ev.candidateTracks_;
  auto& evt_build_extras = ev.candidateTracksExtra_;
  auto& evt_fit_tracks   = ev.fitTracks_;
  auto& evt_fit_extras   = ev.fitTracksExtra_;
  auto& evt_layer_hits   = ev.layerHits_;

  for (auto&& buildtrack : evt_build_tracks)
  {
    const auto& buildextra = evt_build_extras[buildtrack.label()];
    
    // same for fit and build tracks
    evtID_cFR_     = ievt;
    seedID_cFR_    = buildextra.seedID(); 
    mcTrackID_cFR_ = buildextra.mcTrackID();

    // track parameters
    pt_build_cFR_   = buildtrack.pT();
    ept_build_cFR_  = buildtrack.epT();
    phi_build_cFR_  = buildtrack.momPhi();
    ephi_build_cFR_ = buildtrack.emomPhi();
    eta_build_cFR_  = buildtrack.momEta();
    eeta_build_cFR_ = buildtrack.emomEta();
    
    // hit/layer info
    nHits_build_cFR_           = buildtrack.nFoundHits();
    nLayers_build_cFR_         = buildtrack.nUniqueLayers();
    nHitsMatched_build_cFR_    = buildextra.nHitsMatched();
    fracHitsMatched_build_cFR_ = buildextra.fracHitsMatched();
    lastlyr_build_cFR_         = buildtrack.getLastFoundHitLyr();
    
    // hit info
    const Hit& lasthit = evt_layer_hits[buildtrack.getLastFoundHitLyr()][buildtrack.getLastFoundHitIdx()];
    xhit_build_cFR_ = lasthit.x(); 
    yhit_build_cFR_ = lasthit.y(); 
    zhit_build_cFR_ = lasthit.z(); 

    // chi2 info
    hitchi2_build_cFR_ = buildtrack.chi2(); 
    helixchi2_build_cFR_ = buildextra.helixChi2();

    // stored dphi
    dphi_build_cFR_ = buildextra.dPhi();
    
    // cmssw match?
    cmsswID_build_cFR_ = buildextra.cmsswTrackID();
    if (cmsswID_build_cFR_ >= 0) // matched track to cmssw 
    {
      cmsswmask_build_cFR_ = 1; 
    }
    else 
    {
      if (Config::inclusiveShorts) 
      {
	if      (cmsswID_build_cFR_ ==  -1 || cmsswID_build_cFR_ ==  -5 ||
		 cmsswID_build_cFR_ ==  -8 || cmsswID_build_cFR_ ==  -9)  
	{
	  cmsswmask_build_cFR_ = 0;
	}
	else if (cmsswID_build_cFR_ ==  -2 || cmsswID_build_cFR_ == -10 ||
		 cmsswID_build_cFR_ == -11)
	{
	  cmsswmask_build_cFR_ = 2; 
	}
	else // mcID == -3,-4,-6,-7,-12,-13
	{
	  cmsswmask_build_cFR_ = -1;
	}
      }
      else // only count long tracks
      {
	if      (cmsswID_build_cFR_ == -1 || cmsswID_build_cFR_ == -9) 
	{
	  cmsswmask_build_cFR_ = 0;
	}
	else // mcID == -2,-3,-4,-5,-6,-7,-8,-10,-11,-12,-13
	{
	  cmsswmask_build_cFR_ = -1; 
	}
      }
    } // end check over not matched
    
    if (cmsswmask_build_cFR_ == 1) // matched track to cmssw
    {
      const auto& cmsswtrack = evt_cmssw_tracks[cmsswID_build_cFR_];
      const auto& cmsswextra = evt_cmssw_extras[cmsswtrack.label()];

      seedID_cmssw_build_cFR_ = cmsswextra.seedID();

      x_cmssw_build_cFR_ = cmsswtrack.x();
      y_cmssw_build_cFR_ = cmsswtrack.y();
      z_cmssw_build_cFR_ = cmsswtrack.z();
      
      pt_cmssw_build_cFR_  = cmsswtrack.pT(); 
      phi_cmssw_build_cFR_ = cmsswtrack.momPhi();
      eta_cmssw_build_cFR_ = cmsswtrack.momEta();

      nHits_cmssw_build_cFR_   = cmsswtrack.nFoundHits(); 
      nLayers_cmssw_build_cFR_ = cmsswtrack.nUniqueLayers();
      lastlyr_cmssw_build_cFR_ = cmsswtrack.getLastFoundHitLyr();

      // duplicate info
      duplmask_build_cFR_   = buildextra.isDuplicate(); 
      iTkMatches_build_cFR_ = buildextra.duplicateID();
    }
    else // unmatched cmsswtracks ... put -99 for all reco values to denote unmatched
    {
      seedID_cmssw_build_cFR_ = -99;

      x_cmssw_build_cFR_ = -2000;
      y_cmssw_build_cFR_ = -2000;
      z_cmssw_build_cFR_ = -2000;
      
      pt_cmssw_build_cFR_  = -99;
      phi_cmssw_build_cFR_ = -99;
      eta_cmssw_build_cFR_ = -99;

      nHits_cmssw_build_cFR_   = -99;
      nLayers_cmssw_build_cFR_ = -99;
      lastlyr_cmssw_build_cFR_ = -99;

      duplmask_build_cFR_   = -1;
      iTkMatches_build_cFR_ = -99;
    }
  
    // ensure there is a fit track to mess with
    if (buildToFitMap_.count(buildtrack.label()))
    {
      const auto& fittrack = evt_fit_tracks[buildToFitMap_[buildtrack.label()]];
      const auto& fitextra = evt_fit_extras[fittrack.label()];

      // track parameters
      pt_fit_cFR_   = fittrack.pT();
      ept_fit_cFR_  = fittrack.epT();
      phi_fit_cFR_  = fittrack.momPhi();
      ephi_fit_cFR_ = fittrack.emomPhi();
      eta_fit_cFR_  = fittrack.momEta();
      eeta_fit_cFR_ = fittrack.emomEta();
    
      // hit/layer info
      nHits_fit_cFR_           = fittrack.nFoundHits();
      nLayers_fit_cFR_         = fittrack.nUniqueLayers();
      nHitsMatched_fit_cFR_    = fitextra.nHitsMatched();
      fracHitsMatched_fit_cFR_ = fitextra.fracHitsMatched();
      lastlyr_fit_cFR_         = fittrack.getLastFoundHitLyr();
    
      // hit info
      const Hit& lasthit = evt_layer_hits[fittrack.getLastFoundHitLyr()][fittrack.getLastFoundHitIdx()];
      xhit_fit_cFR_ = lasthit.x(); 
      yhit_fit_cFR_ = lasthit.y(); 
      zhit_fit_cFR_ = lasthit.z(); 

      // chi2 info
      hitchi2_fit_cFR_ = fittrack.chi2(); 
      helixchi2_fit_cFR_ = fitextra.helixChi2();

      // stored dphi
      dphi_fit_cFR_ = fitextra.dPhi();
    
      // cmssw match?
      cmsswID_fit_cFR_ = fitextra.cmsswTrackID();
      if (cmsswID_fit_cFR_ >= 0) // matched track to cmssw 
      {
	cmsswmask_fit_cFR_ = 1; 
      }
      else 
      {
	if (Config::inclusiveShorts) 
	{
	  if      (cmsswID_fit_cFR_ ==  -1 || cmsswID_fit_cFR_ ==  -5 ||
		   cmsswID_fit_cFR_ ==  -8 || cmsswID_fit_cFR_ ==  -9)  
	  {
	    cmsswmask_fit_cFR_ = 0;
	  }
	  else if (cmsswID_fit_cFR_ ==  -2 || cmsswID_fit_cFR_ == -10 ||
		   cmsswID_fit_cFR_ == -11)
	  {
	    cmsswmask_fit_cFR_ = 2; 
	  }
	  else // mcID == -3,-4,-6,-7,-12,-13
	  {
	    cmsswmask_fit_cFR_ = -1;
	  }
	}
	else // only count long tracks
	{
	  if      (cmsswID_fit_cFR_ == -1 || cmsswID_fit_cFR_ == -9) 
	  {
	    cmsswmask_fit_cFR_ = 0;
	  }
	  else // mcID == -2,-3,-4,-5,-6,-7,-8,-10,-11,-12,-13
	  {
	    cmsswmask_fit_cFR_ = -1; 
	  }
	}
      } // end check over not matched
      
      if (cmsswmask_fit_cFR_ == 1) // matched track to cmssw
      {
	const auto& cmsswtrack = evt_cmssw_tracks[cmsswID_fit_cFR_];
	const auto& cmsswextra = evt_cmssw_extras[cmsswtrack.label()];
	
	seedID_cmssw_fit_cFR_ = cmsswextra.seedID();
	
	x_cmssw_fit_cFR_ = cmsswtrack.x();
	y_cmssw_fit_cFR_ = cmsswtrack.y();
	z_cmssw_fit_cFR_ = cmsswtrack.z();
      
	pt_cmssw_fit_cFR_  = cmsswtrack.pT(); 
	phi_cmssw_fit_cFR_ = cmsswtrack.momPhi();
	eta_cmssw_fit_cFR_ = cmsswtrack.momEta();
	
	nHits_cmssw_fit_cFR_   = cmsswtrack.nFoundHits(); 
	nLayers_cmssw_fit_cFR_ = cmsswtrack.nUniqueLayers();
	lastlyr_cmssw_fit_cFR_ = cmsswtrack.getLastFoundHitLyr();

	// duplicate info
	duplmask_fit_cFR_   = fitextra.isDuplicate(); 
	iTkMatches_fit_cFR_ = fitextra.duplicateID();
      }
      else // unmatched cmsswtracks ... put -99 for all reco values to denote unmatched
      {
	seedID_cmssw_fit_cFR_ = -99;
	
	x_cmssw_fit_cFR_ = -2000;
	y_cmssw_fit_cFR_ = -2000;
	z_cmssw_fit_cFR_ = -2000;
	
	pt_cmssw_fit_cFR_  = -99;
	phi_cmssw_fit_cFR_ = -99;
	eta_cmssw_fit_cFR_ = -99;

	nHits_cmssw_fit_cFR_   = -99;
	nLayers_cmssw_fit_cFR_ = -99;
	lastlyr_cmssw_fit_cFR_ = -99;
	
	duplmask_fit_cFR_   = -1;
	iTkMatches_fit_cFR_ = -99;
      }
    }
    else // no fit track to match to a build track! 
    {
      pt_fit_cFR_   = -100;
      ept_fit_cFR_  = -100;
      phi_fit_cFR_  = -100;
      ephi_fit_cFR_ = -100;
      eta_fit_cFR_  = -100;
      eeta_fit_cFR_ = -100;
    
      nHits_fit_cFR_           = -100;
      nLayers_fit_cFR_         = -100;
      nHitsMatched_fit_cFR_    = -100;
      fracHitsMatched_fit_cFR_ = -100;
      lastlyr_fit_cFR_         = -100;
    
      xhit_fit_cFR_ = -3000;
      yhit_fit_cFR_ = -3000;
      zhit_fit_cFR_ = -3000;

      hitchi2_fit_cFR_   = -100;
      helixchi2_fit_cFR_ = -100;
      dphi_fit_cFR_      = -100;
    
      cmsswID_fit_cFR_   = -100;
      cmsswmask_fit_cFR_ = -2;
      
      seedID_cmssw_fit_cFR_ = -100;
	
      x_cmssw_fit_cFR_ = -3000;
      y_cmssw_fit_cFR_ = -3000;
      z_cmssw_fit_cFR_ = -3000;
      
      pt_cmssw_fit_cFR_  = -100;
      phi_cmssw_fit_cFR_ = -100;
      eta_cmssw_fit_cFR_ = -100;

      nHits_cmssw_fit_cFR_   = -100;
      nLayers_cmssw_fit_cFR_ = -100;
      lastlyr_cmssw_fit_cFR_ = -100;
      
      duplmask_fit_cFR_   = -2;
      iTkMatches_fit_cFR_ = -100;
    }

    cmsswfrtree_->Fill();
  } 
}

void TTreeValidation::saveTTrees() 
{  
  std::lock_guard<std::mutex> locker(glock_); 
  f_->cd();
  f_->Write();
}             

#endif
