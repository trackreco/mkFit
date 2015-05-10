#include "TTreeValidation.h"
#ifndef NO_ROOT

static bool sortByHitsChi2(const Track& cand1, const Track& cand2)
{
  if (cand1.nHits()==cand2.nHits()) return cand1.chi2()<cand2.chi2();
  return cand1.nHits()>cand2.nHits();
}

static bool tracksByID(const Track& t1, const Track& t2)
{
  return t1.mcTrackID()<t2.mcTrackID();
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
  efftree_->Branch("pt_mc",&pt_mc_);
  efftree_->Branch("pt_seed",&pt_seed_);
  efftree_->Branch("ept_seed",&ept_seed_);
  efftree_->Branch("pt_build",&pt_build_);
  efftree_->Branch("ept_build",&ept_build_);
  efftree_->Branch("pt_fit",&pt_fit_);
  efftree_->Branch("ept_fit",&ept_fit_);

  efftree_->Branch("pz_mc",&pz_mc_);
  efftree_->Branch("pz_seed",&pz_seed_);
  efftree_->Branch("epz_seed",&epz_seed_);
  efftree_->Branch("pz_build",&pz_build_);
  efftree_->Branch("epz_build",&epz_build_);
  efftree_->Branch("pz_fit",&pz_fit_);
  efftree_->Branch("epz_fit",&epz_fit_);

  efftree_->Branch("phi_mc",&phi_mc_);
  efftree_->Branch("phi_seed",&phi_seed_);
  efftree_->Branch("ephi_seed",&ephi_seed_);
  efftree_->Branch("phi_build",&phi_build_);
  efftree_->Branch("ephi_build",&ephi_build_);
  efftree_->Branch("phi_fit",&phi_fit_);
  efftree_->Branch("ephi_fit",&ephi_fit_);

  efftree_->Branch("eta_mc",&eta_mc_);
  efftree_->Branch("eta_seed",&eta_seed_);
  efftree_->Branch("eeta_seed",&eeta_seed_);
  efftree_->Branch("eta_build",&eta_build_);
  efftree_->Branch("eeta_build",&eeta_build_);
  efftree_->Branch("eta_fit",&eta_fit_);
  efftree_->Branch("eeta_fit",&eeta_fit_);

  efftree_->Branch("nHits_mc",&nHits_mc_);
  efftree_->Branch("nHits_seed",&nHits_seed_);
  efftree_->Branch("nHits_build",&nHits_build_);
  efftree_->Branch("nHits_fit",&nHits_fit_);

  efftree_->Branch("chi2_seed",&chi2_seed_);
  efftree_->Branch("chi2_build",&chi2_build_);
  efftree_->Branch("chi2_fit",&chi2_fit_);

  // fake rate validation
  efftree_ = new TTree("efftree","efftree");
  efftree_->Branch("pt_seed",&pt_seed_);
  efftree_->Branch("pt_build",&pt_build_);
  efftree_->Branch("pt_fit",&pt_fit_);

  efftree_->Branch("pz_seed",&pz_seed_);
  efftree_->Branch("pz_build",&pz_build_);
  efftree_->Branch("pz_fit",&pz_fit_);

  efftree_->Branch("phi_seed",&phi_seed_);
  efftree_->Branch("phi_build",&phi_build_);
  efftree_->Branch("phi_fit",&phi_fit_);

  efftree_->Branch("eta_seed",&eta_seed_);
  efftree_->Branch("eta_build",&eta_build_);
  efftree_->Branch("eta_fit",&eta_fit_);

  efftree_->Branch("nHits_seed",&nHits_seed_);
  efftree_->Branch("nHits_build",&nHits_build_);
  efftree_->Branch("nHits_fit",&nHits_fit_);

  efftree_->Branch("chi2_seed",&chi2_seed_);
  efftree_->Branch("chi2_build",&chi2_build_);
  efftree_->Branch("chi2_fit",&chi2_fit_);
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
  mapSimToTks(evt_seed_tracks,simToSeedMap_);
  mapSimToTks(evt_build_tracks,simToBuildMap_);
  mapSimToTks(evt_fit_tracks,simToFitMap_);
}

void TTreeValidation::fillEffTree(const TrackVec& evt_sim_tracks, const TrackVec& evt_seed_tracks, const TrackVec& evt_build_tracks, const TrackVec& evt_fit_tracks){
  std::lock_guard<std::mutex> locker(glock_);

  for (int i = 0; i < evt_sim_tracks.size(); i++){
    pt_mc_    = evt_sim_tracks[i].pt();
    pz_mc_    = evt_sim_tracks[i].pz();
    phi_mc_   = evt_sim_tracks[i].momPhi();
    eta_mc_   = evt_sim_tracks[i].momEta();
    nHits_mc_ = evt_sim_tracks[i].nHits();
    
    if (simToSeedMap_.count(i)){ // recoToSim match : save reco track that first fills map (which by default is the best one --> most hits, lowest chi2)
      pt_seed_    = simToSeedMap_[i][0].pt(); 
      ept_seed_   = simToSeedMap_[i][0].ept();

      pz_seed_    = simToSeedMap_[i][0].pz(); 
      epz_seed_   = simToSeedMap_[i][0].epz();
      
      phi_seed_   = simToSeedMap_[i][0].momPhi(); 
      ephi_seed_  = simToSeedMap_[i][0].emomPhi();
      
      eta_seed_   = simToSeedMap_[i][0].momEta(); 
      eeta_seed_  = simToSeedMap_[i][0].emomEta();

      nHits_seed_ = simToSeedMap_[i][0].nHits();
      chi2_seed_  = simToSeedMap_[i][0].chi2();
    }
    else{ // unmatched simTracks ... put -99 for all reco values to denote unmatched
      pt_seed_    = -99;
      ept_seed_   = -99;

      pz_seed_    = -99;
      epz_seed_   = -99;
     
      phi_seed_   = -99;
      ephi_seed_  = -99;
      
      eta_seed_   = -99;
      eeta_seed_  = -99;
     
      nHits_seed_ = -99;
      chi2_seed_ = -99;
    }

    if (simToBuildMap_.count(i)){ // recoToSim match : save reco track that first fills map (which by default is the best one --> most hits, lowest chi2)
      pt_build_    = simToBuildMap_[i][0].pt(); 
      ept_build_   = simToBuildMap_[i][0].ept();

      pz_build_    = simToBuildMap_[i][0].pz(); 
      epz_build_   = simToBuildMap_[i][0].epz();
      
      phi_build_   = simToBuildMap_[i][0].momPhi(); 
      ephi_build_  = simToBuildMap_[i][0].emomPhi();
      
      eta_build_   = simToBuildMap_[i][0].momEta(); 
      eeta_build_  = simToBuildMap_[i][0].emomEta();

      nHits_build_ = simToBuildMap_[i][0].nHits();
      chi2_build_  = simToBuildMap_[i][0].chi2();
    }
    else{ // unmatched simTracks ... put -99 for all reco values to denote unmatched
      pt_build_    = -99;
      ept_build_   = -99;

      pz_build_    = -99;
      epz_build_   = -99;
     
      phi_build_   = -99;
      ephi_build_  = -99;
      
      eta_build_   = -99;
      eeta_build_  = -99;
     
      nHits_build_ = -99;
      chi2_build_  = -99;
    }

    if (simToFitMap_.count(i)){ // recoToSim match : save reco track that first fills map (which by default is the best one --> most hits, lowest chi2)
      pt_fit_    = simToFitMap_[i][0].pt(); 
      ept_fit_   = simToFitMap_[i][0].ept();

      pz_fit_    = simToFitMap_[i][0].pz(); 
      epz_fit_   = simToFitMap_[i][0].epz();
      
      phi_fit_   = simToFitMap_[i][0].momPhi(); 
      ephi_fit_  = simToFitMap_[i][0].emomPhi();
      
      eta_fit_   = simToFitMap_[i][0].momEta(); 
      eeta_fit_  = simToFitMap_[i][0].emomEta();

      nHits_fit_ = simToFitMap_[i][0].nHits();
      chi2_fit_  = simToFitMap_[i][0].chi2();
    }
    else{ // unmatched simTracks ... put -99 for all reco values to denote unmatched
      pt_fit_    = -99;
      ept_fit_   = -99;

      pz_fit_    = -99;
      epz_fit_   = -99;
     
      phi_fit_   = -99;
      ephi_fit_  = -99;
      
      eta_fit_   = -99;
      eeta_fit_  = -99;
     
      nHits_fit_ = -99;
      chi2_fit_  = -99;
    }
    efftree_->Fill();
  }
}

void TTreeValidation::mapSimToTks(TrackVec& evt_tracks, simToTkMap& simTkMap_){
  //  std::lock_guard<std::mutex> locker(glock_); 
  for (int i = 0; i < evt_tracks.size(); i++){
    evt_tracks[i].setMCTrackID();
  }
  std::sort(evt_tracks.begin(), evt_tracks.end(), tracksByID);

  TrkIter tkIter;
  for (tkIter = evt_tracks.begin(); tkIter != evt_tracks.end(); ++tkIter){
    simTkMap_[(*tkIter).mcTrackID()].push_back((*tkIter));
  }
  
  simToTkMap::iterator tkMapIter;
  for (tkMapIter = simTkMap_.begin(); tkMapIter != simTkMap_.end(); ++tkMapIter){
    if (tkMapIter->second.size() < 2){ // no duplicates, only one matched is by default best one
      continue;
    }
    else{ //sort duplicates to keep the best one --> most hits, lowest chi2
      std::sort(tkMapIter->second.begin(), tkMapIter->second.end(), sortByHitsChi2); 
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
