// N.B. Mask assignments
// --> mcmask_[reco] == 1, "associated" reco to sim track [possible duplmask_[reco] == 1,0] {eff and FR}
// --> mcmask_[reco] == 0, "unassociated" reco to sim track. by definition no duplicates (no reco to associate to sim tracks!) [possible duplmask_[reco] == 2 {eff and FR}]
// --> mcmask_[reco] == -1, "no matching seed to build/fit" track, therefore no build/fit track to match sim! [possible duplmask_[reco] == -1] {FR only} 

// --> nTkMatches_[reco] > 1,   n reco tracks associated to the same sim track ID {eff only}
// --> nTkMatches_[reco] == 1,  1 reco track associated to single sim track ID {eff only}
// --> nTkMatches_[reco] == -99, no reco to sim match {eff only}

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

#include "TTreeValidation.h"
#include "BinInfoUtils.h"
#include "Event.h"
#include "Simulation.h"
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
  // configuration storage
  configtree_ = new TTree("configtree","configtree");
  configtree_->Branch("simtime",&simtime_);
  configtree_->Branch("segtime",&segtime_);
  configtree_->Branch("seedtime",&seedtime_);
  configtree_->Branch("buildtime",&buildtime_);
  configtree_->Branch("fittime",&fittime_);
  configtree_->Branch("hlvtime",&hlvtime_);
  configtree_->Branch("Ntracks",&Ntracks_);
  configtree_->Branch("Nevents",&Nevents_);
  configtree_->Branch("nPhiPart",&nPhiPart_);
  configtree_->Branch("nPhiFactor",&nPhiFactor_);
  configtree_->Branch("nEtaPart",&nEtaPart_);
  configtree_->Branch("etaDet",&etaDet_);
  configtree_->Branch("nlayers_per_seed",&nlayers_per_seed_);

  // build validation
  tree_br_ = new TTree("tree_br","tree_br");
  tree_br_->Branch("layer",&layer_,"layer/i");
  tree_br_->Branch("branches",&branches_,"branches/i");
  tree_br_->Branch("cands",&cands_,"cands/i");
  
  // efficiency validation
  efftree_ = new TTree("efftree","efftree");
  efftree_->Branch("evtID",&evtID_eff_);
  efftree_->Branch("mcID",&mcID_eff_);

  efftree_->Branch("pt_mc_gen",&pt_mc_gen_eff_);
  efftree_->Branch("phi_mc_gen",&phi_mc_gen_eff_);
  efftree_->Branch("eta_mc_gen",&eta_mc_gen_eff_);
  efftree_->Branch("nHits_mc",&nHits_mc_eff_);
  
  efftree_->Branch("x_hit_mc_reco",&x_hit_mc_reco_eff_);
  efftree_->Branch("y_hit_mc_reco",&y_hit_mc_reco_eff_);
  efftree_->Branch("z_hit_mc_reco",&z_hit_mc_reco_eff_);
  
  efftree_->Branch("x_vrx_mc_gen",&x_vrx_mc_gen_eff_);
  efftree_->Branch("y_vrx_mc_gen",&y_vrx_mc_gen_eff_);
  efftree_->Branch("z_vrx_mc_gen",&z_vrx_mc_gen_eff_);

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

void TTreeValidation::fillBuildTree(const unsigned int layer, const unsigned int branches, const unsigned int cands)
{
  std::lock_guard<std::mutex> locker(glock_);
  layer_    = layer;
  branches_ = branches;
  cands_    = cands;
  tree_br_->Fill();
}

void TTreeValidation::makeSimToTksMaps(TrackVec& evt_seed_tracks, TrackVec& evt_build_tracks, TrackVec& evt_fit_tracks){
  // set mcTkIDs... and sort by each (simTracks set in order by default!)
  std::lock_guard<std::mutex> locker(glock_);
  simToSeedMap_.clear();
  simToBuildMap_.clear();
  simToFitMap_.clear();

  mapSimToTks(evt_seed_tracks,simToSeedMap_);
  mapSimToTks(evt_build_tracks,simToBuildMap_);
  mapSimToTks(evt_fit_tracks,simToFitMap_);
}

void TTreeValidation::mapSimToTks(TrackVec& evt_tracks, simToTksMap& simTkMap){
  //  std::lock_guard<std::mutex> locker(glock_); 

  for (auto&& track : evt_tracks){
    track.setMCTrackIDInfo();
    simTkMap[track.mcTrackID()].push_back(&track);
  }

  for (auto&& simTkMatches : simTkMap){
    if (simTkMatches.second.size() < 2) {
      simTkMatches.second[0]->setMCDuplicateInfo(0,bool(false));
    }
    else if (simTkMatches.second[0]->mcTrackID() == 999999){ // do not bother sorting fakes. nobody likes fakes.
      continue;
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

void TTreeValidation::makeSeedToTkMaps(const TrackVec& evt_build_tracks, const TrackVec& evt_fit_tracks){
  std::lock_guard<std::mutex> locker(glock_);
  
  seedToBuildMap_.clear();
  seedToFitMap_.clear();

  mapSeedToTk(evt_build_tracks,seedToBuildMap_);
  mapSeedToTk(evt_fit_tracks,seedToFitMap_);
}

void TTreeValidation::mapSeedToTk(const TrackVec& evt_tracks, seedToTkMap& seedTkMap){
  //  std::lock_guard<std::mutex> locker(glock_); 

  for (auto&& track : evt_tracks){
    seedTkMap[track.seedID()] = &track;
  }
}

void TTreeValidation::fillEffTree(const TrackVec& evt_sim_tracks, const unsigned int& ievt){
  std::lock_guard<std::mutex> locker(glock_);

  for (auto&& simtrack : evt_sim_tracks){
    evtID_eff_ = ievt;
    mcID_eff_  = simtrack.mcTrackID();

    pt_mc_gen_eff_  = simtrack.pt();
    phi_mc_gen_eff_ = simtrack.momPhi();
    eta_mc_gen_eff_ = simtrack.momEta();
    nHits_mc_eff_   = simtrack.nHits();

    x_hit_mc_reco_eff_.clear();
    y_hit_mc_reco_eff_.clear();
    z_hit_mc_reco_eff_.clear();
    for (auto&& simhit : simtrack.hitsVector()){
      x_hit_mc_reco_eff_.push_back(simhit.x());
      y_hit_mc_reco_eff_.push_back(simhit.y());
      z_hit_mc_reco_eff_.push_back(simhit.z());
    }
    
    x_vrx_mc_gen_eff_ = simtrack.x();
    y_vrx_mc_gen_eff_ = simtrack.y();
    z_vrx_mc_gen_eff_ = simtrack.z();

    // matched seed track
    if (simToSeedMap_.count(mcID_eff_)){ // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToSeedMap_[matched SimID][first element in vector]
      mcmask_seed_eff_  = 1; // quick logic for matched

      // use this to access correct sim track layer params
      const float layer = simToSeedMap_[mcID_eff_][0]->hitsVector().back().layer(); // last layer seed ended up on
      const SVector6 & initLayPrms = simtrack.initParamsVector()[layer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps

      pt_mc_seed_eff_  = std::sqrt(getRad2(initLayPrms.At(3),initLayPrms.At(4)));
      phi_mc_seed_eff_ = getPhi(initLayPrms.At(3),initLayPrms.At(4));
      eta_mc_seed_eff_ = getEta(pt_mc_seed_eff_,initLayPrms.At(5));

      pt_seed_eff_   = simToSeedMap_[mcID_eff_][0]->pt(); 
      ept_seed_eff_  = simToSeedMap_[mcID_eff_][0]->ept();
      phi_seed_eff_  = simToSeedMap_[mcID_eff_][0]->momPhi(); 
      ephi_seed_eff_ = simToSeedMap_[mcID_eff_][0]->emomPhi();
      eta_seed_eff_  = simToSeedMap_[mcID_eff_][0]->momEta(); 
      eeta_seed_eff_ = simToSeedMap_[mcID_eff_][0]->emomEta();
   
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
      
      duplmask_seed_eff_ = 2; // mask means unmatched sim track
      nTkMatches_seed_eff_    = -99; // unmatched
    }

    // matched build track
    if (simToBuildMap_.count(mcID_eff_)){ // recoToSim match : save best match --> most hits, lowest chi2, i.e. simToBuildMap_[matched SimID][first element in vector]
      mcmask_build_eff_  = 1; // quick logic for matched

      // use this to access correct sim track layer params
      const float layer = simToBuildMap_[mcID_eff_][0]->hitsVector().back().layer(); // last layer build track ended up on
      const SVector6 & initLayPrms = simtrack.initParamsVector()[layer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps

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
      mcmask_fit_eff_  = 1; // quick logic for matched

      // use this to access correct sim track layer params
      const float layer = simToFitMap_[mcID_eff_][0]->hitsVector().back().layer(); // last layer fit track ended up on
      const SVector6 & initLayPrms = simtrack.initParamsVector()[layer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps

      pt_mc_fit_eff_  = std::sqrt(getRad2(initLayPrms.At(3),initLayPrms.At(4)));
      phi_mc_fit_eff_ = getPhi(initLayPrms.At(3),initLayPrms.At(4));
      eta_mc_fit_eff_ = getEta(pt_mc_fit_eff_,initLayPrms.At(5));

      pt_fit_eff_   = simToFitMap_[mcID_eff_][0]->pt(); 
      ept_fit_eff_  = simToFitMap_[mcID_eff_][0]->ept();
      phi_fit_eff_  = simToFitMap_[mcID_eff_][0]->momPhi(); 
      ephi_fit_eff_ = simToFitMap_[mcID_eff_][0]->emomPhi();
      eta_fit_eff_  = simToFitMap_[mcID_eff_][0]->momEta(); 
      eeta_fit_eff_ = simToFitMap_[mcID_eff_][0]->emomEta();
   
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

void TTreeValidation::fillFakeRateTree(const TrackVec& evt_sim_tracks, const TrackVec& evt_seed_tracks, const unsigned int& ievt){
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
      const SVector6 & initLayPrms = evt_sim_tracks[mcID_seed_FR_].initParamsVector()[layer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps
      
      pt_mc_seed_FR_    = std::sqrt(getRad2(initLayPrms.At(3),initLayPrms.At(4))); // again, pt of mc truth at the layer the seed ends
      phi_mc_seed_FR_   = getPhi(initLayPrms.At(3),initLayPrms.At(4));
      eta_mc_seed_FR_   = getEta(pt_mc_seed_FR_,initLayPrms.At(5));
      nHits_mc_seed_FR_ = evt_sim_tracks[mcID_seed_FR_].nHits();

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
	const SVector6 & initLayPrms = evt_sim_tracks[mcID_build_FR_].initParamsVector()[layer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps
      
	pt_mc_build_FR_    = std::sqrt(getRad2(initLayPrms.At(3),initLayPrms.At(4))); // again, pt of mc truth at the layer the seed ends
	phi_mc_build_FR_   = getPhi(initLayPrms.At(3),initLayPrms.At(4));
	eta_mc_build_FR_   = getEta(pt_mc_build_FR_,initLayPrms.At(5));
	nHits_mc_build_FR_ = evt_sim_tracks[mcID_build_FR_].nHits();

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
	const SVector6 & initLayPrms = evt_sim_tracks[mcID_fit_FR_].initParamsVector()[layer]; //--> can do this as all sim tracks pass through each layer once, and are stored in order... will need to fix this once we have loopers/overlaps
   
	pt_mc_fit_FR_    = std::sqrt(getRad2(initLayPrms.At(3),initLayPrms.At(4))); // again, pt of mc truth at the layer the seed ends
	phi_mc_fit_FR_   = getPhi(initLayPrms.At(3),initLayPrms.At(4));
	eta_mc_fit_FR_   = getEta(pt_mc_fit_FR_,initLayPrms.At(5));
	nHits_mc_fit_FR_ = evt_sim_tracks[mcID_fit_FR_].nHits();

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

void TTreeValidation::fillConfigTree(const unsigned int & ntracks, const unsigned int & nevents, const std::vector<double>& ticks){
  simtime_   = ticks[0];
  segtime_   = ticks[1];
  seedtime_  = ticks[2];
  buildtime_ = ticks[3];
  fittime_   = ticks[4];
  hlvtime_   = ticks[5];

  Ntracks_ = ntracks;
  Nevents_ = nevents;

  nPhiPart_   = Config::nPhiPart;
  nPhiFactor_ = Config::nPhiFactor;
  nEtaPart_   = Config::nEtaPart;
  etaDet_     = Config::etaDet;

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
  
  configtree_->Fill();
}

void TTreeValidation::saveTTrees() {  
  std::lock_guard<std::mutex> locker(glock_); 
  f_->cd();
  f_->Write();
  f_->Close();
}             

#endif
