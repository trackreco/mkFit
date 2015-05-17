// N.B. Mask assignments
// --> mcmask_[reco] == 1, "associated" reco to sim track [possible duplmask_[reco] == 1,0] {eff and FR}
// --> mcmask_[reco] == 0, "unassociated" reco to sim track. by definition no duplicates (no reco to associate to sim tracks!) [possible duplmask_[reco] == 2 {eff and FR}]
// --> mcmask_[reco] == -1, "no matching seed to build/fit" track, therefore no build/fit track to match sim! [possible duplmask_[reco] == -1] {FR only} 

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
#ifndef NO_ROOT

static bool sortByHitsChi2(const Track* cand1, const Track* cand2)
{
  if (cand1->nHits()==cand2->nHits()) return cand1->chi2()<cand2->chi2();
  return cand1->nHits()>cand2->nHits();
}

TTreeValidation::TTreeValidation(std::string fileName)
{
  std::lock_guard<std::mutex> locker(glock_);
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
  
  efftree_->Branch("mcmask_seed",&mcmask_seed_eff_);
  efftree_->Branch("mcmask_build",&mcmask_build_eff_);
  efftree_->Branch("mcmask_fit",&mcmask_fit_eff_);

  efftree_->Branch("pt_mc",&pt_mc_eff_);
  efftree_->Branch("pt_seed",&pt_seed_eff_);
  efftree_->Branch("ept_seed",&ept_seed_eff_);
  efftree_->Branch("pt_build",&pt_build_eff_);
  efftree_->Branch("ept_build",&ept_build_eff_);
  efftree_->Branch("pt_fit",&pt_fit_eff_);
  efftree_->Branch("ept_fit",&ept_fit_eff_);

  efftree_->Branch("pz_mc",&pz_mc_eff_);
  efftree_->Branch("pz_seed",&pz_seed_eff_);
  efftree_->Branch("epz_seed",&epz_seed_eff_);
  efftree_->Branch("pz_build",&pz_build_eff_);
  efftree_->Branch("epz_build",&epz_build_eff_);
  efftree_->Branch("pz_fit",&pz_fit_eff_);
  efftree_->Branch("epz_fit",&epz_fit_eff_);

  efftree_->Branch("phi_mc",&phi_mc_eff_);
  efftree_->Branch("phi_seed",&phi_seed_eff_);
  efftree_->Branch("ephi_seed",&ephi_seed_eff_);
  efftree_->Branch("phi_build",&phi_build_eff_);
  efftree_->Branch("ephi_build",&ephi_build_eff_);
  efftree_->Branch("phi_fit",&phi_fit_eff_);
  efftree_->Branch("ephi_fit",&ephi_fit_eff_);

  efftree_->Branch("eta_mc",&eta_mc_eff_);
  efftree_->Branch("eta_seed",&eta_seed_eff_);
  efftree_->Branch("eeta_seed",&eeta_seed_eff_);
  efftree_->Branch("eta_build",&eta_build_eff_);
  efftree_->Branch("eeta_build",&eeta_build_eff_);
  efftree_->Branch("eta_fit",&eta_fit_eff_);
  efftree_->Branch("eeta_fit",&eeta_fit_eff_);

  efftree_->Branch("chi2_seed",&chi2_seed_eff_);
  efftree_->Branch("chi2_build",&chi2_build_eff_);
  efftree_->Branch("chi2_fit",&chi2_fit_eff_);

  efftree_->Branch("nHits_mc",&nHits_mc_eff_);
  efftree_->Branch("nHits_seed",&nHits_seed_eff_);
  efftree_->Branch("nHits_build",&nHits_build_eff_);
  efftree_->Branch("nHits_fit",&nHits_fit_eff_);

  efftree_->Branch("nHitsMatched_seed",&nHitsMatched_seed_eff_);
  efftree_->Branch("nHitsMatched_build",&nHitsMatched_build_eff_);
  efftree_->Branch("nHitsMatched_fit",&nHitsMatched_fit_eff_);

  efftree_->Branch("duplmask_seed",&duplmask_seed_eff_);
  efftree_->Branch("duplmask_build",&duplmask_build_eff_);
  efftree_->Branch("duplmask_fit",&duplmask_fit_eff_);

  efftree_->Branch("nDup_seed",&nDup_seed_eff_);
  efftree_->Branch("nDup_build",&nDup_build_eff_);
  efftree_->Branch("nDup_fit",&nDup_fit_eff_);

  // fake rate validation
  fakeratetree_ = new TTree("fakeratetree","fakeratetree");

  fakeratetree_->Branch("evtID",&evtID_FR_);
  fakeratetree_->Branch("seedID",&seedID_FR_);

  fakeratetree_->Branch("seedmask_build",&seedmask_build_FR_);
  fakeratetree_->Branch("seedmask_fit",&seedmask_fit_FR_);

  fakeratetree_->Branch("pt_seed",&pt_seed_FR_);
  fakeratetree_->Branch("ept_seed",&ept_seed_FR_);
  fakeratetree_->Branch("pt_build",&pt_build_FR_);
  fakeratetree_->Branch("ept_build",&ept_build_FR_);
  fakeratetree_->Branch("pt_fit",&pt_fit_FR_);
  fakeratetree_->Branch("ept_fit",&ept_fit_FR_);

  fakeratetree_->Branch("pz_seed",&pz_seed_FR_);
  fakeratetree_->Branch("epz_seed",&epz_seed_FR_);
  fakeratetree_->Branch("pz_build",&pz_build_FR_);
  fakeratetree_->Branch("epz_build",&epz_build_FR_);
  fakeratetree_->Branch("pz_fit",&pz_fit_FR_);
  fakeratetree_->Branch("epz_fit",&epz_fit_FR_);

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

  fakeratetree_->Branch("chi2_seed",&chi2_seed_FR_);
  fakeratetree_->Branch("chi2_build",&chi2_build_FR_);
  fakeratetree_->Branch("chi2_fit",&chi2_fit_FR_);

  fakeratetree_->Branch("nHits_seed",&nHits_seed_FR_);
  fakeratetree_->Branch("nHits_build",&nHits_build_FR_);
  fakeratetree_->Branch("nHits_fit",&nHits_fit_FR_);

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

  fakeratetree_->Branch("pz_mc_seed",&pz_mc_seed_FR_);
  fakeratetree_->Branch("pz_mc_build",&pz_mc_build_FR_);
  fakeratetree_->Branch("pz_mc_fit",&pz_mc_fit_FR_);

  fakeratetree_->Branch("phi_mc_seed",&phi_mc_seed_FR_);
  fakeratetree_->Branch("phi_mc_build",&phi_mc_build_FR_);
  fakeratetree_->Branch("phi_mc_fit",&phi_mc_fit_FR_);

  fakeratetree_->Branch("eta_mc_seed",&eta_mc_seed_FR_);
  fakeratetree_->Branch("eta_mc_build",&eta_mc_build_FR_);
  fakeratetree_->Branch("eta_mc_fit",&eta_mc_fit_FR_);

  fakeratetree_->Branch("nHitsMatched_seed",&nHitsMatched_seed_FR_);
  fakeratetree_->Branch("nHitsMatched_build",&nHitsMatched_build_FR_);
  fakeratetree_->Branch("nHitsMatched_fit",&nHitsMatched_fit_FR_);

  fakeratetree_->Branch("nHits_mc_seed",&nHits_mc_seed_FR_);
  fakeratetree_->Branch("nHits_mc_build",&nHits_mc_build_FR_);
  fakeratetree_->Branch("nHits_mc_fit",&nHits_mc_fit_FR_);

  fakeratetree_->Branch("duplmask_seed",&duplmask_seed_FR_);
  fakeratetree_->Branch("duplmask_build",&duplmask_build_FR_);
  fakeratetree_->Branch("duplmask_fit",&duplmask_fit_FR_);

  fakeratetree_->Branch("iDup_seed",&iDup_seed_FR_);
  fakeratetree_->Branch("iDup_build",&iDup_build_FR_);
  fakeratetree_->Branch("iDup_fit",&iDup_fit_FR_);
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

    pt_mc_eff_    = simtrack.pt();
    pz_mc_eff_    = simtrack.pz();
    phi_mc_eff_   = simtrack.momPhi();
    eta_mc_eff_   = simtrack.momEta();

    nHits_mc_eff_ = simtrack.nHits();

    // matched seed track
    if (simToSeedMap_.count(mcID_eff_)){ // recoToSim match : save best match --> most hits, lowest chi2
      mcmask_seed_eff_  = 1; // quick logic for matched

      pt_seed_eff_    = simToSeedMap_[mcID_eff_][0]->pt(); 
      ept_seed_eff_   = simToSeedMap_[mcID_eff_][0]->ept();
      pz_seed_eff_    = simToSeedMap_[mcID_eff_][0]->pz(); 
      epz_seed_eff_   = simToSeedMap_[mcID_eff_][0]->epz();
      phi_seed_eff_   = simToSeedMap_[mcID_eff_][0]->momPhi(); 
      ephi_seed_eff_  = simToSeedMap_[mcID_eff_][0]->emomPhi();
      eta_seed_eff_   = simToSeedMap_[mcID_eff_][0]->momEta(); 
      eeta_seed_eff_  = simToSeedMap_[mcID_eff_][0]->emomEta();

      chi2_seed_eff_         = -10; //simToSeedMap_[mcID_eff_][0]->chi2(); // currently not implemented
      nHits_seed_eff_        = simToSeedMap_[mcID_eff_][0]->nHits();
      nHitsMatched_seed_eff_ = simToSeedMap_[mcID_eff_][0]->nHitsMatched();
      duplmask_seed_eff_     = simToSeedMap_[mcID_eff_][0]->isDuplicate(); 
      nDup_seed_eff_         = simToSeedMap_[mcID_eff_].size(); // n reco matches to this sim track.
    }
    else{ // unmatched simTracks ... put -99 for all reco values to denote unmatched
      mcmask_seed_eff_  = 0; // quick logic for not matched

      pt_seed_eff_    = -99;
      ept_seed_eff_   = -99;
      pz_seed_eff_    = -99;
      epz_seed_eff_   = -99;
      phi_seed_eff_   = -99;
      ephi_seed_eff_  = -99;
      eta_seed_eff_   = -99;
      eeta_seed_eff_  = -99;

      chi2_seed_eff_         = -99;
      nHits_seed_eff_        = -99;
      nHitsMatched_seed_eff_ = -99;
      duplmask_seed_eff_     = 2; // mask means unmatched sim track therefore should enter denom8
      nDup_seed_eff_         = -99; // unmatched
    }

    // matched build track
    if (simToBuildMap_.count(mcID_eff_)){ 
      mcmask_build_eff_  = 1; // quick logic for matched

      pt_build_eff_    = simToBuildMap_[mcID_eff_][0]->pt(); 
      ept_build_eff_   = simToBuildMap_[mcID_eff_][0]->ept();
      pz_build_eff_    = simToBuildMap_[mcID_eff_][0]->pz(); 
      epz_build_eff_   = simToBuildMap_[mcID_eff_][0]->epz();
      phi_build_eff_   = simToBuildMap_[mcID_eff_][0]->momPhi(); 
      ephi_build_eff_  = simToBuildMap_[mcID_eff_][0]->emomPhi();
      eta_build_eff_   = simToBuildMap_[mcID_eff_][0]->momEta(); 
      eeta_build_eff_  = simToBuildMap_[mcID_eff_][0]->emomEta();

      chi2_build_eff_         = -10; //simToBuildMap_[mcID_eff_][0]->chi2(); // currently not implemented
      nHits_build_eff_        = simToBuildMap_[mcID_eff_][0]->nHits();
      nHitsMatched_build_eff_ = simToBuildMap_[mcID_eff_][0]->nHitsMatched();
      duplmask_build_eff_     = simToBuildMap_[mcID_eff_][0]->isDuplicate(); 
      nDup_build_eff_         = simToBuildMap_[mcID_eff_].size(); // n reco matches to this sim track.
    }
    else{ // unmatched simTracks ... put -99 for all reco values to denote unmatched
      mcmask_build_eff_  = 0; // quick logic for not matched

      pt_build_eff_    = -99;
      ept_build_eff_   = -99;
      pz_build_eff_    = -99;
      epz_build_eff_   = -99;
      phi_build_eff_   = -99;
      ephi_build_eff_  = -99;
      eta_build_eff_   = -99;
      eeta_build_eff_  = -99;

      chi2_build_eff_         = -99;
      nHits_build_eff_        = -99;
      nHitsMatched_build_eff_ = -99;
      duplmask_build_eff_     = 2;
      nDup_build_eff_         = -99; // unmatched
    }

    // matched fit track
    if (simToFitMap_.count(mcID_eff_)){ // recoToSim match : save best match --> most hits, lowest chi2
      mcmask_fit_eff_  = 1; // quick logic for matched

      pt_fit_eff_    = simToFitMap_[mcID_eff_][0]->pt(); 
      ept_fit_eff_   = simToFitMap_[mcID_eff_][0]->ept();
      pz_fit_eff_    = simToFitMap_[mcID_eff_][0]->pz(); 
      epz_fit_eff_   = simToFitMap_[mcID_eff_][0]->epz();
      phi_fit_eff_   = simToFitMap_[mcID_eff_][0]->momPhi(); 
      ephi_fit_eff_  = simToFitMap_[mcID_eff_][0]->emomPhi();
      eta_fit_eff_   = simToFitMap_[mcID_eff_][0]->momEta(); 
      eeta_fit_eff_  = simToFitMap_[mcID_eff_][0]->emomEta();

      chi2_fit_eff_         = -10; //simToFitMap_[mcID_eff_][0]->chi2(); // currently not implemented
      nHits_fit_eff_        = simToFitMap_[mcID_eff_][0]->nHits();
      nHitsMatched_fit_eff_ = simToFitMap_[mcID_eff_][0]->nHitsMatched();
      duplmask_fit_eff_     = simToFitMap_[mcID_eff_][0]->isDuplicate(); 
      nDup_fit_eff_         = simToFitMap_[mcID_eff_].size(); // n reco matches to this sim track.
    }
    else{ // unmatched simTracks ... put -99 for all reco values to denote unmatched
      mcmask_fit_eff_  = 0; // quick logic for not matched

      pt_fit_eff_    = -99;
      ept_fit_eff_   = -99;
      pz_fit_eff_    = -99;
      epz_fit_eff_   = -99;
      phi_fit_eff_   = -99;
      ephi_fit_eff_  = -99;
      eta_fit_eff_   = -99;
      eeta_fit_eff_  = -99;

      chi2_fit_eff_         = -99;
      nHits_fit_eff_        = -99;
      nHitsMatched_fit_eff_ = -99;
      duplmask_fit_eff_     = 2;
      nDup_fit_eff_         = -99; // unmatched
    }

    efftree_->Fill(); // fill it once per sim track!
  }
}

void TTreeValidation::fillFakeRateTree(const TrackVec& evt_sim_tracks, const TrackVec& evt_seed_tracks, const unsigned int& ievt){
  std::lock_guard<std::mutex> locker(glock_);
  
  // these operations are expensive, do this once
  std::vector<float> simpt;
  std::vector<float> simpz;
  std::vector<float> simphi;
  std::vector<float> simeta;
  std::vector<unsigned int> simnHits;

  for (auto&& simtrack : evt_sim_tracks){ // assume no sorting on simtracks! if so, then resort sim tracks by ID
    simpt.push_back(simtrack.pt());
    simpz.push_back(simtrack.pz());
    simphi.push_back(simtrack.momPhi());
    simeta.push_back(simtrack.momEta());
    simnHits.push_back(simtrack.nHits());
  }

  for (auto&& seedtrack : evt_seed_tracks){
    evtID_FR_  = ievt;
    seedID_FR_ = seedtrack.seedID();

    // seed info
    pt_seed_FR_   = seedtrack.pt();
    ept_seed_FR_  = seedtrack.ept();
    pz_seed_FR_   = seedtrack.pz();
    epz_seed_FR_  = seedtrack.epz();
    phi_seed_FR_  = seedtrack.momPhi();
    ephi_seed_FR_ = seedtrack.emomPhi();
    eta_seed_FR_  = seedtrack.momEta();
    eeta_seed_FR_ = seedtrack.emomEta();

    chi2_seed_FR_  = -10; // seedtrack.chi2(); --> not yet implemented
    nHits_seed_FR_ = seedtrack.nHits();
    
    // sim info for seed track
    mcID_seed_FR_ = seedtrack.mcTrackID();
    if (mcID_seed_FR_ != 999999){
      mcmask_seed_FR_ = 1; // matched track to sim
      
      pt_mc_seed_FR_  = simpt[mcID_seed_FR_];
      pz_mc_seed_FR_  = simpz[mcID_seed_FR_];
      phi_mc_seed_FR_ = simphi[mcID_seed_FR_];
      eta_mc_seed_FR_ = simeta[mcID_seed_FR_];

      nHitsMatched_seed_FR_ = seedtrack.nHitsMatched();
      nHits_mc_seed_FR_ = simnHits[mcID_seed_FR_];

      duplmask_seed_FR_ = seedtrack.isDuplicate();
      iDup_seed_FR_     = seedtrack.duplicateID(); // ith duplicate seed track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"      
    }
    else{
      mcmask_seed_FR_ = 0;   // fake track (unmatched track)
          
      // -99 for all sim info for reco tracks not associated to reco tracks
      pt_mc_seed_FR_  = -99;
      pz_mc_seed_FR_  = -99;
      phi_mc_seed_FR_ = -99;
      eta_mc_seed_FR_ = -99;

      nHitsMatched_seed_FR_ = -99;
      nHits_mc_seed_FR_ = -99;

      duplmask_seed_FR_ = 2; // see notation above      
      iDup_seed_FR_     = -99;  
    }

    //==========================//
    
    // fill build information if track still alive
    if (seedToBuildMap_.count(seedID_FR_)){
      seedmask_build_FR_ = 1; // quick logic

      pt_build_FR_   = seedToBuildMap_[seedID_FR_]->pt();
      ept_build_FR_  = seedToBuildMap_[seedID_FR_]->ept();
      pz_build_FR_   = seedToBuildMap_[seedID_FR_]->pz();
      epz_build_FR_  = seedToBuildMap_[seedID_FR_]->epz();
      phi_build_FR_  = seedToBuildMap_[seedID_FR_]->momPhi();
      ephi_build_FR_ = seedToBuildMap_[seedID_FR_]->emomPhi();
      eta_build_FR_  = seedToBuildMap_[seedID_FR_]->momEta();
      eeta_build_FR_ = seedToBuildMap_[seedID_FR_]->emomEta();
      
      chi2_build_FR_  = seedToBuildMap_[seedID_FR_]->chi2();
      nHits_build_FR_ = seedToBuildMap_[seedID_FR_]->nHits();

      mcID_build_FR_  = seedToBuildMap_[seedID_FR_]->mcTrackID();
      if (mcID_build_FR_ != 999999){ // build track matched to seed and sim 
	mcmask_build_FR_ = 1; // matched track to sim
	
	pt_mc_build_FR_  = simpt[mcID_build_FR_];
	pz_mc_build_FR_  = simpz[mcID_build_FR_];
	phi_mc_build_FR_ = simphi[mcID_build_FR_];
	eta_mc_build_FR_ = simeta[mcID_build_FR_];
	
	nHits_mc_build_FR_     = simnHits[mcID_build_FR_];

	nHitsMatched_build_FR_ = seedToBuildMap_[seedID_FR_]->nHitsMatched();
	duplmask_build_FR_ = seedToBuildMap_[seedID_FR_]->isDuplicate();
	iDup_build_FR_     = seedToBuildMap_[seedID_FR_]->duplicateID(); // ith duplicate build track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"
      }
      else{ // build track matched only to seed not to sim
	mcmask_build_FR_ = 0;   // fake track (unmatched track)
	
	// -99 for all sim info for reco tracks not associated to reco tracks
	pt_mc_build_FR_  = -99;
	pz_mc_build_FR_  = -99;
	phi_mc_build_FR_ = -99;
	eta_mc_build_FR_ = -99;
	
	nHits_mc_build_FR_     = -99;

	nHitsMatched_build_FR_ = -99;
	duplmask_build_FR_ = 2;
	iDup_build_FR_     = -99;
      } // matched seed to build, not build to sim
    }

    else { // seed has no matching build track (therefore no matching sim to build track)
      seedmask_build_FR_ = 0; // quick logic

      // -100 for all reco info as no actual build track for this seed
      pt_build_FR_   = -100;
      ept_build_FR_  = -100;
      pz_build_FR_   = -100;
      epz_build_FR_  = -100;
      phi_build_FR_  = -100;
      ephi_build_FR_ = -100;
      eta_build_FR_  = -100;
      eeta_build_FR_ = -100;
      
      chi2_build_FR_  = -100; 
      nHits_build_FR_ = -100;

      // keep -100 for all sim variables as no such reco exists for this seed
      mcmask_build_FR_ = -1; // do not want to count towards build FR
      mcID_build_FR_   = -100;
	
      pt_mc_build_FR_  = -100;
      pz_mc_build_FR_  = -100;
      phi_mc_build_FR_ = -100;
      eta_mc_build_FR_ = -100;
      
      nHits_mc_build_FR_ = -100;

      nHitsMatched_build_FR_ = -100;
      duplmask_build_FR_ = -1;
      iDup_build_FR_     = -100;
    }

    //============================// fit tracks

    if (seedToFitMap_.count(seedID_FR_)){
      seedmask_fit_FR_ = 1; // quick logic

      pt_fit_FR_   = seedToFitMap_[seedID_FR_]->pt();
      ept_fit_FR_  = seedToFitMap_[seedID_FR_]->ept();
      pz_fit_FR_   = seedToFitMap_[seedID_FR_]->pz();
      epz_fit_FR_  = seedToFitMap_[seedID_FR_]->epz();
      phi_fit_FR_  = seedToFitMap_[seedID_FR_]->momPhi();
      ephi_fit_FR_ = seedToFitMap_[seedID_FR_]->emomPhi();
      eta_fit_FR_  = seedToFitMap_[seedID_FR_]->momEta();
      eeta_fit_FR_ = seedToFitMap_[seedID_FR_]->emomEta();
      
      chi2_fit_FR_  = -10; // seedToFitMap_[seedID_FR_]->chi2(); --> not yet implemented
      nHits_fit_FR_ = seedToFitMap_[seedID_FR_]->nHits();

      mcID_fit_FR_  = seedToFitMap_[seedID_FR_]->mcTrackID();
      if (mcID_fit_FR_ != 999999){ // fit track matched to seed and sim 
	mcmask_fit_FR_ = 1; // matched track to sim
	
	pt_mc_fit_FR_  = simpt[mcID_fit_FR_];
	pz_mc_fit_FR_  = simpz[mcID_fit_FR_];
	phi_mc_fit_FR_ = simphi[mcID_fit_FR_];
	eta_mc_fit_FR_ = simeta[mcID_fit_FR_];
	
	nHits_mc_fit_FR_     = simnHits[mcID_fit_FR_];

	nHitsMatched_fit_FR_ = seedToFitMap_[seedID_FR_]->nHitsMatched();
	duplmask_fit_FR_ = seedToFitMap_[seedID_FR_]->isDuplicate();
	iDup_fit_FR_     = seedToFitMap_[seedID_FR_]->duplicateID(); // ith duplicate fit track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"
      }
      else{ // fit track matched only to seed not to sim
	mcmask_fit_FR_ = 0;   // fake track (unmatched track)
	
	// -99 for all sim info for reco tracks not associated to reco tracks
	pt_mc_fit_FR_  = -99;
	pz_mc_fit_FR_  = -99;
	phi_mc_fit_FR_ = -99;
	eta_mc_fit_FR_ = -99;
	
	nHits_mc_fit_FR_     = -99;

	nHitsMatched_fit_FR_ = -99;
	duplmask_fit_FR_ = 2;
	iDup_fit_FR_     = -99;
      } // matched seed to fit, not fit to sim
    }

    else { // seed has no matching fit track (therefore no matching sim to fit track)
      seedmask_fit_FR_ = 0; // quick logic

      // -100 for all reco info as no actual fit track for this seed
      pt_fit_FR_   = -100;
      ept_fit_FR_  = -100;
      pz_fit_FR_   = -100;
      epz_fit_FR_  = -100;
      phi_fit_FR_  = -100;
      ephi_fit_FR_ = -100;
      eta_fit_FR_  = -100;
      eeta_fit_FR_ = -100;

      chi2_fit_FR_  = -100; 
      nHits_fit_FR_ = -100;

      // keep -100 for all sim variables as no such reco exists for this seed
      mcmask_fit_FR_ = -1; // do not want to count towards fit FR
      mcID_fit_FR_   = -100;
	
      pt_mc_fit_FR_  = -100;
      pz_mc_fit_FR_  = -100;
      phi_mc_fit_FR_ = -100;
      eta_mc_fit_FR_ = -100;
      
      nHits_mc_fit_FR_ = -100;
      
      nHitsMatched_fit_FR_ = -100;
      duplmask_fit_FR_ = -1;
      iDup_fit_FR_     = -100;
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
 
  configtree_->Fill();
}

void TTreeValidation::saveTTrees() {  
  std::lock_guard<std::mutex> locker(glock_); 
  f_->cd();
  f_->Write();
  f_->Close();
}             

#endif
