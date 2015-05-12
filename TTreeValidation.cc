#include "TTreeValidation.h"
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

  // fake rate validation
  fakeseedtree_ = new TTree("fakeseedtree","fakeseedtree");
  fakeseedtree_->Branch("mask_seed_fake",&mask_seed_fake_);
  fakeseedtree_->Branch("mask_seed_duplicate",&mask_seed_duplicate_);
  fakeseedtree_->Branch("pt_seed",&pt_seed_fake_);
  fakeseedtree_->Branch("pz_seed",&pz_seed_fake_);
  fakeseedtree_->Branch("phi_seed",&phi_seed_fake_);
  fakeseedtree_->Branch("eta_seed",&eta_seed_fake_);
  fakeseedtree_->Branch("nHits_seed",&nHits_seed_fake_);
  fakeseedtree_->Branch("chi2_seed",&chi2_seed_fake_);
  fakeseedtree_->Branch("evt_seed",&evt_seed_fake_);

  fakebuildtree_ = new TTree("fakebuildtree","fakebuildtree");
  fakebuildtree_->Branch("mask_build_fake",&mask_build_fake_);
  fakebuildtree_->Branch("mask_build_duplicate",&mask_build_duplicate_);
  fakebuildtree_->Branch("pt_build",&pt_build_fake_);
  fakebuildtree_->Branch("pz_build",&pz_build_fake_);
  fakebuildtree_->Branch("phi_build",&phi_build_fake_);
  fakebuildtree_->Branch("eta_build",&eta_build_fake_);
  fakebuildtree_->Branch("nHits_build",&nHits_build_fake_);
  fakebuildtree_->Branch("chi2_build",&chi2_build_fake_);
  fakebuildtree_->Branch("evt_build",&evt_build_fake_);

  fakefittree_ = new TTree("fakefittree","fakefittree");
  fakefittree_->Branch("mask_fit_fake",&mask_fit_fake_);
  fakefittree_->Branch("mask_fit_duplicate",&mask_fit_duplicate_);
  fakefittree_->Branch("pt_fit",&pt_fit_fake_);
  fakefittree_->Branch("pz_fit",&pz_fit_fake_);
  fakefittree_->Branch("phi_fit",&phi_fit_fake_);
  fakefittree_->Branch("eta_fit",&eta_fit_fake_);
  fakefittree_->Branch("nHits_fit",&nHits_fit_fake_);
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

void TTreeValidation::fillEffTree(const TrackVec& evt_sim_tracks, const unsigned int& ievt){
  std::lock_guard<std::mutex> locker(glock_);

  for (int i = 0; i < evt_sim_tracks.size(); i++){
    pt_mc_eff_    = evt_sim_tracks[i].pt();
    pz_mc_eff_    = evt_sim_tracks[i].pz();
    phi_mc_eff_   = evt_sim_tracks[i].momPhi();
    eta_mc_eff_   = evt_sim_tracks[i].momEta();
    nHits_mc_eff_ = evt_sim_tracks[i].nHits();
    evt_mc_eff_   = ievt;
  
    if (simToSeedMap_.count(i)){ // recoToSim match : save best match --> most hits, lowest chi2
      pt_seed_eff_    =  simToSeedMap_[i][0]->pt(); 
      ept_seed_eff_   =  simToSeedMap_[i][0]->ept();
      pz_seed_eff_    =  simToSeedMap_[i][0]->pz(); 
      epz_seed_eff_   =  simToSeedMap_[i][0]->epz();
      phi_seed_eff_   =  simToSeedMap_[i][0]->momPhi(); 
      ephi_seed_eff_  =  simToSeedMap_[i][0]->emomPhi();
      eta_seed_eff_   =  simToSeedMap_[i][0]->momEta(); 
      eeta_seed_eff_  =  simToSeedMap_[i][0]->emomEta();
      nHits_seed_eff_ =  simToSeedMap_[i][0]->nHits();
      nHitsMatched_seed_eff_ =  simToSeedMap_[i][0]->nHitsMatched();
      chi2_seed_eff_  =  simToSeedMap_[i][0]->chi2();
      evt_seed_eff_   = ievt;
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
      evt_seed_eff_   = -99;
    }

    if (simToBuildMap_.count(i)){ // recoToSim match : save best match --> most hits, lowest chi2
      pt_build_eff_    =  simToBuildMap_[i][0]->pt(); 
      ept_build_eff_   =  simToBuildMap_[i][0]->ept();
      pz_build_eff_    =  simToBuildMap_[i][0]->pz(); 
      epz_build_eff_   =  simToBuildMap_[i][0]->epz();
      phi_build_eff_   =  simToBuildMap_[i][0]->momPhi(); 
      ephi_build_eff_  =  simToBuildMap_[i][0]->emomPhi();
      eta_build_eff_   =  simToBuildMap_[i][0]->momEta(); 
      eeta_build_eff_  =  simToBuildMap_[i][0]->emomEta();
      nHits_build_eff_ =  simToBuildMap_[i][0]->nHits();
      nHitsMatched_build_eff_ =  simToBuildMap_[i][0]->nHitsMatched();
      chi2_build_eff_  =  simToBuildMap_[i][0]->chi2();
      evt_build_eff_   = ievt;
    }
    else{ // unmatched simTracks ... put -99 for all reco values to denote unmatched
      pt_build_eff_    = -99;
      ept_build_eff_   = -99;
      pz_build_eff_    = -99;
      epz_build_eff_   = -99;
      phi_build_eff_   = -99;
      ephi_build_eff_  = -99;
      eta_build_eff_   = -99;
      eeta_build_eff_  = -99;
      nHits_build_eff_ = -99;
      nHitsMatched_build_eff_ = -99;
      chi2_build_eff_  = -99;
      evt_build_eff_   = -99;
    }

    if (simToFitMap_.count(i)){ // recoToSim match : save best match --> most hits, lowest chi2
      pt_fit_eff_    =  simToFitMap_[i][0]->pt(); 
      ept_fit_eff_   =  simToFitMap_[i][0]->ept();
      pz_fit_eff_    =  simToFitMap_[i][0]->pz(); 
      epz_fit_eff_   =  simToFitMap_[i][0]->epz();
      phi_fit_eff_   =  simToFitMap_[i][0]->momPhi(); 
      ephi_fit_eff_  =  simToFitMap_[i][0]->emomPhi();
      eta_fit_eff_   =  simToFitMap_[i][0]->momEta(); 
      eeta_fit_eff_  =  simToFitMap_[i][0]->emomEta();
      nHits_fit_eff_ =  simToFitMap_[i][0]->nHits();
      nHitsMatched_fit_eff_ =  simToFitMap_[i][0]->nHitsMatched();
      chi2_fit_eff_  =  simToFitMap_[i][0]->chi2();
      evt_fit_eff_   = ievt;
    }
    else{ // unmatched simTracks ... put -99 for all reco values to denote unmatched
      pt_fit_eff_    = -99;
      ept_fit_eff_   = -99;
      pz_fit_eff_    = -99;
      epz_fit_eff_   = -99;
      phi_fit_eff_   = -99;
      ephi_fit_eff_  = -99;
      eta_fit_eff_   = -99;
      eeta_fit_eff_  = -99;
      nHits_fit_eff_ = -99;
      nHitsMatched_fit_eff_ = -99;
      chi2_fit_eff_  = -99;
      evt_fit_eff_   = -99;
    }

    efftree_->Fill(); // fill it once per sim track!
  }
}

void TTreeValidation::fillFakeTrees(const unsigned int& ievt){
  //  std::lock_guard<std::mutex> locker(glock_);

  fillFakeSeedTree(ievt);
  fillFakeBuildTree(ievt);
  fillFakeFitTree(ievt);
}

void TTreeValidation::fillFakeSeedTree(const unsigned int& ievt){
  std::lock_guard<std::mutex> locker(glock_);

  for (auto&& recPairs : simToSeedMap_){
    unsigned int i = 0; // i is dummy counter ... want to keep track of where you are in vector... 1st element is ref to "best track"
    for (auto&& track : recPairs.second){
      if (track->mcTrackID() != 999999){
	mask_seed_fake_ = 1;
	if (i == 0){
	  mask_seed_duplicate_ = 1; // 1 means "best candidate"
       	}
	else{
	  mask_seed_duplicate_ = 0; // 0 means "not best candidate, but still real reco"
	}
      }
      else{
	mask_seed_fake_ = 0;
	mask_seed_duplicate_ = -1; // means "don't count towards duplicate_ info" 
      }

      pt_seed_fake_    = track->pt();
      pz_seed_fake_    = track->pz();
      phi_seed_fake_   = track->momPhi();
      eta_seed_fake_   = track->momEta();
      nHits_seed_fake_ = track->nHits();
      chi2_seed_fake_  = track->chi2();
      evt_seed_fake_   = ievt;

      i++; // dummy counter
      fakeseedtree_->Fill();
    }
  }
}

void TTreeValidation::fillFakeBuildTree(const unsigned int& ievt){
  std::lock_guard<std::mutex> locker(glock_);
  
  for (auto&& recPairs : simToBuildMap_){
    unsigned int i = 0; // i is dummy counter ... want to keep track of where you are in vector... 1st element is ref to "best track"
    for (auto&& track : recPairs.second){
      if (track->mcTrackID() != 999999){
	mask_build_fake_ = 1;
	if (i == 0){
	  mask_build_duplicate_ = 1; // 1 means "best candidate"
       	}
	else{
	  mask_build_duplicate_ = 0; // 0 means "not best candidate, but still real reco"
	}
      }
      else{
	mask_build_fake_ = 0;
	mask_build_duplicate_ = -1; // means "don't count towards duplicate_ info" 
      }

      pt_build_fake_    = track->pt();
      pz_build_fake_    = track->pz();
      phi_build_fake_   = track->momPhi();
      eta_build_fake_   = track->momEta();
      nHits_build_fake_ = track->nHits();
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
	mask_fit_fake_ = 1;
	if (i == 0){
	  mask_fit_duplicate_ = 1; // 1 means "best candidate"
       	}
	else{
	  mask_fit_duplicate_ = 0; // 0 means "not best candidate, but still real reco"
	}
      }
      else{
	mask_fit_fake_ = 0;
	mask_fit_duplicate_ = -1; // means "don't count towards duplicate_ info" 
      }

      pt_fit_fake_    = track->pt();
      pz_fit_fake_    = track->pz();
      phi_fit_fake_   = track->momPhi();
      eta_fit_fake_   = track->momEta();
      nHits_fit_fake_ = track->nHits();
      chi2_fit_fake_  = track->chi2();
      evt_fit_fake_   = ievt;
      
      i++; // dummy counter
      fakefittree_->Fill();
    }
  }
}

void TTreeValidation::saveTTrees() {  
  std::lock_guard<std::mutex> locker(glock_); 
  f_->cd();
  f_->Write();
  f_->Close();
}             

#endif
