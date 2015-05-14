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

  efftree_->Branch("nHits_mc",&nHits_mc_eff_);
  efftree_->Branch("nHits_seed",&nHits_seed_eff_);
  efftree_->Branch("nHits_build",&nHits_build_eff_);
  efftree_->Branch("nHits_fit",&nHits_fit_eff_);

  efftree_->Branch("chi2_seed",&chi2_seed_eff_);
  efftree_->Branch("chi2_build",&chi2_build_eff_);
  efftree_->Branch("chi2_fit",&chi2_fit_eff_);

  efftree_->Branch("nHitsMatched_seed",&nHitsMatched_seed_eff_);
  efftree_->Branch("nHitsMatched_build",&nHitsMatched_build_eff_);
  efftree_->Branch("nHitsMatched_fit",&nHitsMatched_fit_eff_);

  efftree_->Branch("evt_mc",&evt_mc_eff_);
  efftree_->Branch("evt_seed",&evt_seed_eff_);
  efftree_->Branch("evt_build",&evt_build_eff_);
  efftree_->Branch("evt_fit",&evt_fit_eff_);
  
  efftree_->Branch("nDup_seed",&nDup_seed_eff_);
  efftree_->Branch("nDup_build",&nDup_build_eff_);
  efftree_->Branch("nDup_seed",&nDup_fit_eff_);

  efftree_->Branch("mask_seed",&mask_seed_eff_);
  efftree_->Branch("mask_build",&mask_build_eff_);
  efftree_->Branch("mask_fit",&mask_fit_eff_);

  // fake rate validation
  fakeseedtree_ = new TTree("fakeseedtree","fakeseedtree");
  fakeseedtree_->Branch("mask_seed_real",&mask_seed_real_);
  fakeseedtree_->Branch("mask_seed_duplicate",&mask_seed_duplicate_);
  fakeseedtree_->Branch("pt_seed",&pt_seed_fake_);
  fakeseedtree_->Branch("pz_seed",&pz_seed_fake_);
  fakeseedtree_->Branch("phi_seed",&phi_seed_fake_);
  fakeseedtree_->Branch("eta_seed",&eta_seed_fake_);
  fakeseedtree_->Branch("nHits_seed",&nHits_seed_fake_);
  fakeseedtree_->Branch("nHitsMatched_seed",&nHitsMatched_seed_fake_);
  fakeseedtree_->Branch("chi2_seed",&chi2_seed_fake_);
  fakeseedtree_->Branch("evt_seed",&evt_seed_fake_);

  fakebuildtree_ = new TTree("fakebuildtree","fakebuildtree");
  fakebuildtree_->Branch("mask_build_real",&mask_build_real_);
  fakebuildtree_->Branch("mask_build_duplicate",&mask_build_duplicate_);
  fakebuildtree_->Branch("pt_build",&pt_build_fake_);
  fakebuildtree_->Branch("pz_build",&pz_build_fake_);
  fakebuildtree_->Branch("phi_build",&phi_build_fake_);
  fakebuildtree_->Branch("eta_build",&eta_build_fake_);
  fakebuildtree_->Branch("nHits_build",&nHits_build_fake_);
  fakebuildtree_->Branch("nHitsMatched_build",&nHitsMatched_build_fake_);
  fakebuildtree_->Branch("chi2_build",&chi2_build_fake_);
  fakebuildtree_->Branch("evt_build",&evt_build_fake_);

  fakefittree_ = new TTree("fakefittree","fakefittree");
  fakefittree_->Branch("mask_fit_real",&mask_fit_real_);
  fakefittree_->Branch("mask_fit_duplicate",&mask_fit_duplicate_);
  fakefittree_->Branch("pt_fit",&pt_fit_fake_);
  fakefittree_->Branch("pz_fit",&pz_fit_fake_);
  fakefittree_->Branch("phi_fit",&phi_fit_fake_);
  fakefittree_->Branch("eta_fit",&eta_fit_fake_);
  fakefittree_->Branch("nHits_fit",&nHits_fit_fake_);
  fakefittree_->Branch("nHitsMatched_fit",&nHitsMatched_fit_fake_);
  fakefittree_->Branch("chi2_fit",&chi2_fit_fake_);
  fakefittree_->Branch("evt_fit",&evt_fit_fake_);
}

void TTreeValidation::fillBuildTree(const unsigned int layer, const unsigned int branches, const unsigned int cands)
{
  std::lock_guard<std::mutex> locker(glock_);
  layer_    = layer;
  branches_ = branches;
  cands_    = cands;
  tree_br_->Fill();
}

void TTreeValidation::makeSimToTkMaps(TrackVec& evt_seed_tracks, TrackVec& evt_build_tracks, TrackVec& evt_fit_tracks){
  // set mcTkIDs... and sort by each (simTracks set in order by default!)
  std::lock_guard<std::mutex> locker(glock_);
  simToSeedMap_.clear();
  simToBuildMap_.clear();
  simToFitMap_.clear();

  mapSimToTks(evt_seed_tracks,simToSeedMap_);
  mapSimToTks(evt_build_tracks,simToBuildMap_);
  mapSimToTks(evt_fit_tracks,simToFitMap_);
}

void TTreeValidation::makeSeedToTkMaps(TrackVec& evt_seed_tracks, TrackVec& evt_build_tracks, TrackVec& evt_fit_tracks){
  std::lock_guard<std::mutex> locker(glock_);
  
  seedToSeedMap_.clear();
  seedToBuildMap_.clear();
  seedToFitMap_.clear();

  mapSeedToTks(evt_seed_tracks,seedToSeedMap_);
  mapSeedToTks(evt_build_tracks,seedToBuildMap_);
  mapSeedToTks(evt_fit_tracks,seedToFitMap_);
}

void TTreeValidation::mapSimToTks(TrackVec& evt_tracks, simToTkMap& simTkMap_){
  //  std::lock_guard<std::mutex> locker(glock_); 

  for (auto&& track : evt_tracks){
    track.setMCTrackIDInfo();
    simTkMap_[track.mcTrackID()].push_back(&track);
  }

  for (auto&& simTkMatches : simTkMap_){
    if (simTkMatches.second.size() < 2) {
      continue;
    }
    else if (simTkMatches.second[0]->mcTrackID() == 999999){ // do not bother sorting fakes. nobody likes fakes.
      continue;
    }
    else{ // sort duplicates (ghosts) to keep best one --> most hits, lowest chi2
      std::sort(simTkMatches.second.begin(), simTkMatches.second.end(), sortByHitsChi2); 
    }
  }
}

void TTreeValidation::mapSeedToTks(TrackVec& evt_tracks, simToTkMap& simTkMap_){
  //  std::lock_guard<std::mutex> locker(glock_); 

  for (auto&& track : evt_tracks){
    simTkMap_[track.seedID()].push_back(&track);
  }

  for (auto&& simTkMatches : simTkMap_){
    if (simTkMatches.second.size() < 2) {
      continue;
    }
    else{ // sort duplicates (ghosts) to keep best one --> most hits, lowest chi2 --> could have more than one track per seed!
      std::sort(simTkMatches.second.begin(), simTkMatches.second.end(), sortByHitsChi2); 
    }
  }
}

void TTreeValidation::fillEffTree(const TrackVec& evt_sim_tracks, const unsigned int& ievt){
  std::lock_guard<std::mutex> locker(glock_);

  for (int i = 0; i < evt_sim_tracks.size(); i++){
    /*    x_reco_mc_.clear();  x_gen_mc_.clear();
    y_reco_mc_.clear();  y_gen_mc_.clear();
    z_reco_mc_.clear();  z_gen_mc_.clear();

    for (auto&& hit : evt_sim_tracks[i].hitsVector()){
      x_reco_mc_.push_back(hit.position()[0]);
      y_reco_mc_.push_back(hit.position()[1]);
      z_reco_mc_.push_back(hit.position()[2]);
    }
    for (auto&& hit : evt_sim_tracks[i].initHitsVector()){
      x_gen_mc_.push_back(hit.position()[0]);
      y_gen_mc_.push_back(hit.position()[1]);
      z_gen_mc_.push_back(hit.position()[2]);
    }
    */    

    pt_mc_eff_    = evt_sim_tracks[i].pt();
    pz_mc_eff_    = evt_sim_tracks[i].pz();
    phi_mc_eff_   = evt_sim_tracks[i].momPhi();
    eta_mc_eff_   = evt_sim_tracks[i].momEta();

    nHits_mc_eff_ = evt_sim_tracks[i].nHits();

    mcID_mc_eff_  = i; // i is simTrackID
    evt_mc_eff_   = ievt;
  
    if (simToSeedMap_.count(i)){ // recoToSim match : save best match --> most hits, lowest chi2
      pt_seed_eff_    = simToSeedMap_[i][0]->pt(); 
      ept_seed_eff_   = simToSeedMap_[i][0]->ept();
      pz_seed_eff_    = simToSeedMap_[i][0]->pz(); 
      epz_seed_eff_   = simToSeedMap_[i][0]->epz();
      phi_seed_eff_   = simToSeedMap_[i][0]->momPhi(); 
      ephi_seed_eff_  = simToSeedMap_[i][0]->emomPhi();
      eta_seed_eff_   = simToSeedMap_[i][0]->momEta(); 
      eeta_seed_eff_  = simToSeedMap_[i][0]->emomEta();

      nHits_seed_eff_ = simToSeedMap_[i][0]->nHits();
      nHitsMatched_seed_eff_ = simToSeedMap_[i][0]->nHitsMatched();
      chi2_seed_eff_  = -10; //simToSeedMap_[i][0]->chi2(); // currently not implemented

      mask_seed_eff_  = 1; // quick logic for matched
      nDup_seed_eff_  = simToSeedMap_[i].size(); // n reco matches to this sim track, just no reco info
    }
    else{ // unmatched simTracks ... put -99 for all reco values to denote unmatched
      pt_seed_eff_    = -99;
      ept_seed_eff_   = -99;
      pz_seed_eff_    = -99;
      epz_seed_eff_   = -99;
      phi_seed_eff_   = -99;
      ephi_seed_eff_  = -99;
      eta_seed_eff_   = -99;
      eeta_seed_eff_  = -99;

      nHits_seed_eff_ = -99;
      nHitsMatched_seed_eff_ = -99;
      chi2_seed_eff_  = -99;

      // for all unmatched simTracks, store info on that sim track

      mask_seed_eff_  = 0; // quick logic for not matched
      nDup_seed_eff_  = -1; // unmatched
    }

    efftree_->Fill(); // fill it once per sim track!
  }
}

void TTreeValidation::fillFakeRateTree(const unsigned int& ievt){
  std::lock_guard<std::mutex> locker(glock_);
  
  std::vector<float> simpt;
  std::vector<float> simpz;
  std::vector<float> simphi;
  std::vector<float> simeta;
  std::vector<unsigned int> simnHits;

  for (auto&& simtrack : evt_sim_tracks){
    simpt.push_back(simtrack.pt());
    simpz.push_back(simtrack.pz());
    simphi.push_back(simtrack.phi());
    simeta.push_back(simtrack.eta());
    simnHits.push_back(simtrack.nHits());
  }

  for (auto&& recPairs : seedToSeedMap_){
    unsigned int i = 0; // i is dummy counter ... want to keep track of where you are in vector... 1st element is ref to "best track"
    for (auto&& track : recPairs.second){
      pt_seed_FR_   = track->pt();
      ept_seed_FR_  = track->ept();
      pz_seed_FR_   = track->pz();
      epz_seed_FR_  = track->epz();
      phi_seed_FR_  = track->momPhi();
      ephi_seed_FR_ = track->emomPhi();
      eta_seed_FR_  = track->momEta();
      eeta_seed_FR_ = track->emomEta();

      nHits_seed_FR_ = track->nHits();
      nHitsMatched_seed_FR_ = track->nHitsMatched();
      chi2_seed_FR_  = track->chi2();

      mcID_seed_FR_ = track->mcTrackID();
      if (track->mcTrackID() != 999999){
	mask_seed_FR_ = 1; // matched track
	iDup_seed_FR_ = i; // ith duplicate track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"

	pt_mc_seed_FR_  = simpt[mcID_seed_FR_];
	pz_mc_seed_FR_  = simpz[mcID_seed_FR_];
	phi_mc_seed_FR_ = simphi[mcID_seed_FR_];
	eta_mc_seed_FR_ = simeta[mcID_seed_FR_];
	
	nHits_mc_seed_FR_ = simnHits[mcID_seed_FR_];
      }
      else{
	mask_seed_FR_ = 0;   // fake track (unmatched track)
 	iDup_seed_FR_ = -1;  // means "don't count towards duplicate_ info" 

	// -99 for all sim info for reco tracks not associated to reco tracks
	pt_mc_seed_FR_  = -99;
	pz_mc_seed_FR_  = -99;
	phi_mc_seed_FR_ = -99;
	eta_mc_seed_FR_ = -99;
	
	nHits_mc_seed_FR_ = -99;
      }

      evt_seed_FR_   = ievt;
      fakeseedtree_->Fill();
      i++; // dummy counter
    }
  } // end of seed to seed loop

  for (auto&& recPairs : seedToBuildMap_){
    unsigned int i = 0; // i is dummy counter ... want to keep track of where you are in vector... 1st element is ref to "best track"
    for (auto&& track : recPairs.second){
      pt_build_FR_   = track->pt();
      ept_build_FR_  = track->ept();
      pz_build_FR_   = track->pz();
      epz_build_FR_  = track->epz();
      phi_build_FR_  = track->momPhi();
      ephi_build_FR_ = track->emomPhi();
      eta_build_FR_  = track->momEta();
      eeta_build_FR_ = track->emomEta();

      nHits_build_FR_ = track->nHits();
      nHitsMatched_build_FR_ = track->nHitsMatched();
      chi2_build_FR_  = track->chi2();

      mcID_build_FR_ = track->mcTrackID();
      if (track->mcTrackID() != 999999){
	mask_build_FR_ = 1; // matched track
	iDup_build_FR_ = i; // ith duplicate track, i = 0 "best" match, i > 0 "still matched, real reco, not as good as i-1 track"

	pt_mc_build_FR_  = simpt[mcID_build_FR_];
	pz_mc_build_FR_  = simpz[mcID_build_FR_];
	phi_mc_build_FR_ = simphi[mcID_build_FR_];
	eta_mc_build_FR_ = simeta[mcID_build_FR_];
	
	nHits_mc_build_FR_ = simnHits[mcID_build_FR_];
      }
      else{
	mask_build_FR_ = 0;   // fake track (unmatched track)
 	iDup_build_FR_ = -1;  // means "don't count towards duplicate_ info" 

	// -99 for all sim info for reco tracks not associated to reco tracks
	pt_mc_build_FR_  = -99;
	pz_mc_build_FR_  = -99;
	phi_mc_build_FR_ = -99;
	eta_mc_build_FR_ = -99;
	
	nHits_mc_build_FR_ = -99;
      }

      evt_build_FR_   = ievt;
      fakebuildtree_->Fill();
      i++; // dummy counter
    }
  }
}


void TTreeValidation::fillFakeSeedTree(const unsigned int& ievt){
  std::lock_guard<std::mutex> locker(glock_);

}

void TTreeValidation::fillFakeBuildTree(const unsigned int& ievt){
  std::lock_guard<std::mutex> locker(glock_);
  
  for (auto&& recPairs : simToBuildMap_){
    unsigned int i = 0; // i is dummy counter ... want to keep track of where you are in vector... 1st element is ref to "best track"
    for (auto&& track : recPairs.second){
      if (track->mcTrackID() != 999999){
	mask_build_real_      = 1;
	mask_build_duplicate_ = i; // 0 means "best candidate", i > 0 means "ith duplpicate that is real reco but not as good as ith-1 track"
      }
      else{
	mask_build_real_ = 0;
	mask_build_duplicate_ = -1; // means "don't count towards duplicate_ info" 
      }

      pt_build_fake_    = track->pt();
      pz_build_fake_    = track->pz();
      phi_build_fake_   = track->momPhi();
      eta_build_fake_   = track->momEta();
      nHits_build_fake_ = track->nHits();
      nHitsMatched_build_fake_ = track->nHitsMatched();
      chi2_build_fake_  = track->chi2();
      evt_build_fake_   = ievt;
      
      i++; // dummy counter
      fakebuildtree_->Fill();
    }
  }
}

void TTreeValidation::fillFakeFitTree(const unsigned int& ievt){
  std::lock_guard<std::mutex> locker(glock_);
  
  for (auto&& recPairs : simToFitMap_){
    unsigned int i = 0; // i is dummy counter ... want to keep track of where you are in vector... 1st element is ref to "best track"
    for (auto&& track : recPairs.second){
      if (track->mcTrackID() != 999999){
	mask_fit_real_      = 1;
	mask_fit_duplicate_ = i;
      }
      else{
	mask_fit_real_      = 0;
	mask_fit_duplicate_ = -1; // means "don't count towards duplicate_ info" 
      }

      pt_fit_fake_    = track->pt();
      pz_fit_fake_    = track->pz();
      phi_fit_fake_   = track->momPhi();
      eta_fit_fake_   = track->momEta();
      nHits_fit_fake_ = track->nHits();
      nHitsMatched_fit_fake_ = track->nHitsMatched();
      chi2_fit_fake_  = track->chi2();
      evt_fit_fake_   = ievt;
      
      i++; // dummy counter
      fakefittree_->Fill();
    }
  }
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
