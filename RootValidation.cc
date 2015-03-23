#include "RootValidation.h"
#ifndef NO_ROOT

float getPt(float px, float py) { return sqrt(px*px + py*py); }
float deltaPhi(float phi1, float phi2) {
  float dphi = std::abs(phi1 - phi2);
  if (dphi > TMath::Pi()) { dphi = (2.*TMath::Pi()) - dphi; }
  return dphi;
}
float deltaEta(float eta1, float eta2) { return (eta1 - eta2); }
float deltaR(float phi1, float eta1, float phi2, float eta2) { 
  return sqrt( deltaPhi(phi1,phi2)*deltaPhi(phi1,phi2) + deltaEta(eta1,eta2)*deltaEta(eta1,eta2) );
}

RootValidation::RootValidation(std::string fileName, bool saveTree)
  : savetree_(saveTree)
{
  std::lock_guard<std::mutex> locker(glock_);
  setupHists();

  if (savetree_) {
    f_ = TFile::Open(fileName.c_str(), "recreate");

    // build validation
    buildtree_ = new TTree("buildtree","buildtree");
    buildtree_->Branch("nhits",&tk_nhits_,"nhits/i");
    buildtree_->Branch("chi2",&tk_chi2_,"chi2/F");
    tree_br_ = new TTree("tree_br","tree_br");
    tree_br_->Branch("layer",&layer_,"layer/i");
    tree_br_->Branch("branches",&branches_,"branches/i");
    tree_br_->Branch("cands",&cands_,"cands/i");

    // fit validation
    fittree_ = new TTree("fittree","fittree");
    fittree_->Branch("pt_mc",&pt_mc,"pt_mc");
    fittree_->Branch("pt_fit",&pt_fit,"pt_fit");
    fittree_->Branch("pt_err",&pt_err,"pt_err");

    fittree_->Branch("simHit0_x",&simHit0_x,"simHit0_x");
    fittree_->Branch("simHit0_y",&simHit0_y,"simHit0_y");
    fittree_->Branch("simHit0_z",&simHit0_z,"simHit0_z");
    fittree_->Branch("simHit0_px",&simHit0_px,"simHit0_px");
    fittree_->Branch("simHit0_py",&simHit0_py,"simHit0_py");
    fittree_->Branch("simHit0_pz",&simHit0_pz,"simHit0_pz");
    fittree_->Branch("cfitHit0_x",&cfitHit0_x,"cfitHit0_x");
    fittree_->Branch("cfitHit0_y",&cfitHit0_y,"cfitHit0_y");
    fittree_->Branch("cfitHit0_z",&cfitHit0_z,"cfitHit0_z");
    fittree_->Branch("cfitHit0_px",&cfitHit0_px,"cfitHit0_px");
    fittree_->Branch("cfitHit0_py",&cfitHit0_py,"cfitHit0_py");
    fittree_->Branch("cfitHit0_pz",&cfitHit0_pz,"cfitHit0_pz");
    fittree_->Branch("cfitHit0_xe",&cfitHit0_xe,"cfitHit0_xe");
    fittree_->Branch("cfitHit0_ye",&cfitHit0_ye,"cfitHit0_ye");
    fittree_->Branch("cfitHit0_ze",&cfitHit0_ze,"cfitHit0_ze");
    fittree_->Branch("cfitHit0_pxe",&cfitHit0_pxe,"cfitHit0_pxe");
    fittree_->Branch("cfitHit0_pye",&cfitHit0_pye,"cfitHit0_pye");
    fittree_->Branch("cfitHit0_pze",&cfitHit0_pze,"cfitHit0_pze");

    posTree_ = new TTree("posTree","posTree");
    posTree_->Branch("x_init",&x_init,"x_init");
    posTree_->Branch("x_mc",&x_mc,"x_mc");
    posTree_->Branch("x_mcerr",&x_mcerr,"x_mcerr");
    posTree_->Branch("x_prop",&x_prop,"x_prop");
    posTree_->Branch("x_perr",&x_perr,"x_perr");
    posTree_->Branch("x_update",&x_update,"x_update");
    posTree_->Branch("x_uerr",&x_uerr,"x_uerr");

    posTree_->Branch("y_init",&y_init,"y_init");
    posTree_->Branch("y_mc",&y_mc,"y_mc");
    posTree_->Branch("y_mcerr",&y_mcerr,"y_mcerr");
    posTree_->Branch("y_prop",&y_prop,"y_prop");
    posTree_->Branch("y_perr",&y_perr,"y_perr");
    posTree_->Branch("y_update",&y_update,"y_update");
    posTree_->Branch("y_uerr",&y_uerr,"y_uerr");

    posTree_->Branch("z_init",&z_init,"z_init");
    posTree_->Branch("z_mc",&z_mc,"z_mc");
    posTree_->Branch("z_mcerr",&z_mcerr,"z_mcerr");
    posTree_->Branch("z_prop",&z_prop,"z_prop");
    posTree_->Branch("z_perr",&z_perr,"z_perr");
    posTree_->Branch("z_update",&z_update,"z_update");
    posTree_->Branch("z_uerr",&z_uerr,"z_uerr");

    posTree_->Branch("xy_mcerr",&xy_mcerr,"xy_mcerr");

    posTree_->Branch("r_init",&r_init,"r_init");
    posTree_->Branch("r_mc",&r_mc,"r_mc");
    posTree_->Branch("r_prop",&r_prop,"r_prop");
    posTree_->Branch("r_update",&r_update,"r_update");

    posTree_->Branch("phi_init",&phi_init,"phi_init");
    posTree_->Branch("phi_mc",&phi_mc,"phi_mc");
    posTree_->Branch("phi_mcerr",&phi_mcerr,"phi_mcerr");
    posTree_->Branch("phi_prop",&phi_prop,"phi_prop");
    posTree_->Branch("phi_perr",&phi_perr,"phi_perr");
    posTree_->Branch("phi_update",&phi_update,"phi_update");
    posTree_->Branch("phi_uerr",&phi_uerr,"phi_uerr");

  }
}

void RootValidation::fillSimHists(const TrackVec& evt_sim_tracks)
{
  // these are expensive, only do once per track
  std::vector<float> pt;
  std::vector<float> phi;
  std::vector<float> eta;
  for (auto&& simtrack : evt_sim_tracks) {
    float tkpt = getPt(simtrack.momentum()[0],simtrack.momentum()[1]);
    float tkphi = getPhi(simtrack.momentum()[0], simtrack.momentum()[1]);
    pt.push_back(tkpt);
    phi.push_back(tkphi);
    float tkptrad = sqrt(simtrack.momentum()[0]*simtrack.momentum()[0]+simtrack.momentum()[1]*simtrack.momentum()[1]);
    float tketa = getEta(tkptrad, simtrack.momentum()[2]);
    eta.push_back(tketa);

    HitVec simhits = simtrack.hitsVector();

    for (auto&& simhit : simhits){
      validation_hists2_["detectorzr"]->Fill(simhit.position()[2],simhit.r());
      validation_hists2_["detectorxy"]->Fill(simhit.position()[0],simhit.position()[1]);
    }
  }

  std::lock_guard<std::mutex> locker(glock_);
  for(unsigned int isim_track = 0; isim_track < evt_sim_tracks.size(); ++isim_track){
    validation_hists_["gen_trk_Pt"]->Fill( pt[isim_track] );
    validation_hists_["gen_trk_Px"]->Fill( evt_sim_tracks[isim_track].momentum()[0] );
    validation_hists_["gen_trk_Py"]->Fill( evt_sim_tracks[isim_track].momentum()[1] ); 
    validation_hists_["gen_trk_Pz"]->Fill( evt_sim_tracks[isim_track].momentum()[2] ); 
    validation_hists_["gen_trk_phi"]->Fill( phi[isim_track] );
    validation_hists_["gen_trk_eta"]->Fill( eta[isim_track] );
    
    const HitVec& hits = evt_sim_tracks[isim_track].hitsVector();
    for (auto&& hit : hits){
      float rad = sqrt(hit.position()[0]*hit.position()[0] + hit.position()[1]*hit.position()[1]);
      validation_hists_["gen_hits_rad"]->Fill( rad );
      // Fill histo for layer 3
      if ( (rad > 11.0) && (rad < 13.0) ) {
        validation_hists_["gen_hits_rad_lay3"]->Fill( rad );
      }
      validation_hists_["gen_hits_cov00"]->Fill( hit.error()[0][0] );
      validation_hists_["gen_hits_cov11"]->Fill( hit.error()[1][1] );
    }
    
    float mindR = 999999;
    float mindPhi = 999999;
    const float phii=phi[isim_track];
    const float etai=eta[isim_track];

    auto&& gen_trk_dR(validation_hists_["gen_trk_dR"]);
    auto&& gen_trk_dPhi(validation_hists_["gen_trk_dPhi"]);

    // doubly nested loop over # tracks, do as little as possible in the inner loop
    for( unsigned int jsim_track = 0; jsim_track < evt_sim_tracks.size(); ++jsim_track ){
      if(jsim_track != isim_track) {
        const float phij=phi[jsim_track];
        const float etaj=eta[jsim_track];
        const float drij = deltaR(phii, etai, phij, etaj);
        const float dphiij = deltaPhi(phii, phij);
        mindR=std::min(mindR, drij);
        mindPhi=std::min(mindPhi, dphiij);
        if(jsim_track > isim_track){
          gen_trk_dR->Fill( drij );
          gen_trk_dPhi->Fill( dphiij );
        }
      }
    } 
    validation_hists_["gen_trk_mindR"]->Fill( mindR );
    validation_hists_["gen_trk_mindPhi"]->Fill( mindPhi );
    validation_hists_["gen_trk_nHits"]->Fill( evt_sim_tracks[isim_track].nHits());
  }
  
}

void RootValidation::fillSeedHists(std::vector<HitVec>& seed_pairs, std::vector<HitVec>& seed_triplets, TrackVec& evt_sim_tracks, HitVec& missed_phis, HitVec & missed_etas, std::vector<float> & deta, std::vector<float> & chi2fit, TrackVec & conformalTracks, TrackVec & seedTracks, TrackVec & seedTracksMC){

  // Do these once per track -- sim
  std::vector<float> simpt;
  std::vector<float> simphi;
  std::vector<float> simeta;
  std::vector<float> simtheta;
  for (auto&& track : evt_sim_tracks) {
    simpt.push_back(getPt(track.momentum()[0], track.momentum()[1]));
    simphi.push_back(getPhi(track.momentum()[0], track.momentum()[1]));
    float rad = sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1]);
    simeta.push_back(getEta(rad, track.momentum()[2]));
    simtheta.push_back(atan2(rad,track.momentum()[2]));
  }

  for (auto&& seed_pair : seed_pairs) {
    unsigned int simtrack = seed_pair[1].mcTrackID();
    if(seed_pair[0].mcTrackID() == simtrack){
      // for efficiency studies
      validation_hists_["matchedSeedPair_SimPt"]->Fill(simpt[simtrack]);
      validation_hists_["matchedSeedPair_SimPhi"]->Fill(simphi[simtrack]);
      validation_hists_["matchedSeedPair_SimEta"]->Fill(simeta[simtrack]);
    } // end criterion if block
    else{
      // just count the number of extra seeds at this point
      validation_hists_["extraSeedPair_SimPt"]->Fill(simpt[simtrack]);
      validation_hists_["extraSeedPair_SimPhi"]->Fill(simphi[simtrack]);
      validation_hists_["extraSeedPair_SimEta"]->Fill(simeta[simtrack]);
    }
  }

  for (unsigned int i = 0; i < seed_triplets.size(); i++){
    unsigned int simtrack = seed_triplets[i][1].mcTrackID();
    if( (seed_triplets[i][0].mcTrackID() == simtrack) && (seed_triplets[i][2].mcTrackID() == simtrack) ){
      validation_hists_["matchedSeedTriplet_chi2fit"]->Fill(chi2fit[i]);
    }
    else{
      validation_hists_["extraSeedTriplet_chi2fit"]->Fill(chi2fit[i]);
    }
  }

  for (auto&& seed_triplet : seed_triplets) {
    unsigned int simtrack = seed_triplet[1].mcTrackID();
    if( (seed_triplet[0].mcTrackID() == simtrack) && (seed_triplet[2].mcTrackID() == simtrack) ){
      // for efficiency studies
      validation_hists_["matchedSeedTriplet_SimPt"]->Fill(simpt[simtrack]);
      validation_hists_["matchedSeedTriplet_SimPhi"]->Fill(simphi[simtrack]);
      validation_hists_["matchedSeedTriplet_SimEta"]->Fill(simeta[simtrack]);
    } // end criterion if block
    else{
      // just count the number of extra seeds at this point
      validation_hists_["extraSeedTriplet_SimPt"]->Fill(simpt[simtrack]);
      validation_hists_["extraSeedTriplet_SimPhi"]->Fill(simphi[simtrack]);
      validation_hists_["extraSeedTriplet_SimEta"]->Fill(simeta[simtrack]);
    }
  }

  for (auto&& conformalTrack : conformalTracks){
    unsigned int simtrack = conformalTrack.hitsVector()[1].mcTrackID();
    if( (conformalTrack.hitsVector()[0].mcTrackID() == simtrack) && (conformalTrack.hitsVector()[2].mcTrackID() == simtrack) ){
      float pt  = getPt(conformalTrack.momentum()[0],conformalTrack.momentum()[1]);
      float phi = getPhi(conformalTrack.momentum()[0],conformalTrack.momentum()[1]);
      float theta = atan2(pt,conformalTrack.momentum()[2]);
      float pt_err = sqrt( 
			  conformalTrack.errors()[3][3]*conformalTrack.momentum()[0]*conformalTrack.momentum()[0] +
			  conformalTrack.errors()[4][4]*conformalTrack.momentum()[1]*conformalTrack.momentum()[1] +
			  2*conformalTrack.errors()[3][4]*conformalTrack.momentum()[0]*conformalTrack.momentum()[1] )/pt;
      validation_hists_["conf_pt_pull"]->Fill((pt-simpt[simtrack])/pt_err);
      //      validation_hists_["conf_pt_res"]->Fill((pt-simpt[simtrack])/simpt[simtrack]);
      validation_hists_["conf_phi_res"]->Fill((phi-simphi[simtrack]));
      validation_hists_["conf_invpt_res"]->Fill((1./pt-1./simpt[simtrack]));
      validation_hists_["conf_theta_res"]->Fill(theta-simtheta[simtrack]);
    }
  }

  for (auto&& seedTrack : seedTracks){
    unsigned int simtrack = seedTrack.hitsVector()[1].mcTrackID();
    if( (seedTrack.hitsVector()[0].mcTrackID() == simtrack) && (seedTrack.hitsVector()[2].mcTrackID() == simtrack) ){
      float pt = getPt(seedTrack.momentum()[0],seedTrack.momentum()[1]);
      float pt_err = sqrt( 
			  seedTrack.errors()[3][3]*seedTrack.momentum()[0]*seedTrack.momentum()[0] +
			  seedTrack.errors()[4][4]*seedTrack.momentum()[1]*seedTrack.momentum()[1] +
			  2*seedTrack.errors()[3][4]*seedTrack.momentum()[0]*seedTrack.momentum()[1] )/pt;
      validation_hists_["seed_pt_pull"]->Fill((pt-simpt[simtrack])/pt_err);
      validation_hists_["seed_pt_res"]->Fill((pt-simpt[simtrack]));
      //validation_hists_["seed_pt_res"]->Fill((pt-simpt[simtrack])/simpt[simtrack]);

      validation_hists_["matchedSeed_SimPt"]->Fill(simpt[simtrack]);
      validation_hists_["matchedSeed_SimPhi"]->Fill(simphi[simtrack]);
      validation_hists_["matchedSeed_SimEta"]->Fill(simeta[simtrack]);
    } // end criterion if block

    else{
      // just count the number of extra seeds at this point
      validation_hists_["extraSeed_SimPt"]->Fill(simpt[simtrack]);
      validation_hists_["extraSeed_SimPhi"]->Fill(simphi[simtrack]);
      validation_hists_["extraSeed_SimEta"]->Fill(simeta[simtrack]);
    }
  }

  for (auto&& seedTrackMC : seedTracksMC){
    unsigned int simtrack = seedTrackMC.hitsVector()[1].mcTrackID();
    if( (seedTrackMC.hitsVector()[0].mcTrackID() == simtrack) && (seedTrackMC.hitsVector()[2].mcTrackID() == simtrack) ){
      float pt = getPt(seedTrackMC.momentum()[0],seedTrackMC.momentum()[1]);
      float pt_err = sqrt( 
			  seedTrackMC.errors()[3][3]*seedTrackMC.momentum()[0]*seedTrackMC.momentum()[0] +
			  seedTrackMC.errors()[4][4]*seedTrackMC.momentum()[1]*seedTrackMC.momentum()[1] +
			  2*seedTrackMC.errors()[3][4]*seedTrackMC.momentum()[0]*seedTrackMC.momentum()[1] )/pt;
      validation_hists_["seedMC_pt_pull"]->Fill((pt-simpt[simtrack])/pt_err);
      //      validation_hists_["seedMC_pt_res"]->Fill((pt-simpt[simtrack])/simpt[simtrack]);
      validation_hists_["seedMC_pt_res"]->Fill(1./pt-1./simpt[simtrack]);
    }
  }

  for (auto&& missed_phi : missed_phis){
    unsigned int simtrack = missed_phi.mcTrackID();
    validation_hists_["missed_phi_SimPt"]->Fill(simpt[simtrack]);
    validation_hists_["missed_phi_SimPhi"]->Fill(simphi[simtrack]);
    validation_hists_["missed_phi_SimEta"]->Fill(simeta[simtrack]);
  }

  for (auto&& missed_eta : missed_etas){
    unsigned int simtrack = missed_eta.mcTrackID();
    validation_hists_["missed_eta_SimPt"]->Fill(simpt[simtrack]);
    validation_hists_["missed_eta_SimPhi"]->Fill(simphi[simtrack]);
    validation_hists_["missed_eta_SimEta"]->Fill(simeta[simtrack]);
  }

  for (auto&& d_eta : deta){
    validation_hists_["deta"]->Fill(d_eta);
  }



}

void RootValidation::fillCandidateHists(const TrackVec& evt_track_candidates)
{
  std::lock_guard<std::mutex> locker(glock_);
  //dump candidates
  for (auto&& tkcand : evt_track_candidates) {
    //    if (tkcand.nHits() >= 10){
      validation_hists_["rec_trk_nHits"]->Fill(tkcand.nHits());
      validation_hists_["rec_trk_chi2"]->Fill(tkcand.chi2());
      validation_hists_["rec_trk_phi"]->Fill( getPhi(tkcand.momentum()[0], tkcand.momentum()[1]) );
      validation_hists_["rec_trk_Pt"]->Fill( getPt(tkcand.momentum()[0], tkcand.momentum()[1]) );
      float rad = sqrt(tkcand.momentum()[0]*tkcand.momentum()[0]+tkcand.momentum()[1]*tkcand.momentum()[1]);
      validation_hists_["rec_trk_eta"]->Fill( getEta(rad, tkcand.momentum()[2]) );
      //    }
    if (savetree_) {
      tk_nhits_ = tkcand.nHits();
      tk_chi2_ = tkcand.chi2();
      buildtree_->Fill();
    }
  }
}

void RootValidation::fillAssociationHists(const TrackVec& evt_track_candidates, const TrackVec& evt_sim_tracks){
  std::lock_guard<std::mutex> locker(glock_);
  //setup for assocation; these are dense in simIndex, so use a vector
  std::vector<unsigned int> associated_indices_found_RD(evt_sim_tracks.size()); 
  std::vector<unsigned int> associated_indices_found_SD(evt_sim_tracks.size()); 

  // Do these once per track -- sim
  std::vector<float> simpt;
  std::vector<float> simphi;
  std::vector<float> simeta;
  for (auto&& track : evt_sim_tracks) {
    simpt.push_back(getPt(track.momentum()[0], track.momentum()[1]));
    simphi.push_back(getPhi(track.momentum()[0], track.momentum()[1]));
    float rad = sqrt(track.momentum()[0]*track.momentum()[0]+track.momentum()[1]*track.momentum()[1]);
    simeta.push_back(getEta(rad, track.momentum()[2]));
  }

  for (auto&& tkcand : evt_track_candidates) {

    // get sim Index --> matching only simTrackID, not hitIDs... probably not a big deal for loopers/overlap
    SimTkIDInfo candSimIDInfo = tkcand.SimTrackIDInfo();
    unsigned int simtrack     = candSimIDInfo.first;
    unsigned int nHitsMatched = candSimIDInfo.second;

    float candpt  = getPt(tkcand.momentum()[0],tkcand.momentum()[1]);
    float candphi = getPhi(tkcand.momentum()[0],tkcand.momentum()[1]);
    float candrad = sqrt(tkcand.momentum()[0]*tkcand.momentum()[0]+tkcand.momentum()[1]*tkcand.momentum()[1]);
    float candeta = getEta(candrad,tkcand.momentum()[2]);

    // Check to see if reco track has enough hits matched -- RD
    unsigned int denom_nHits_RD = tkcand.nHits();
    //    if (denom_nHits_RD >= 10){
    if (4*nHitsMatched >= 3*denom_nHits_RD){ // if association criterion is passed, save the info
      if (associated_indices_found_RD[simtrack] == 0){ // currently unmatched simtrack, count it towards efficiency 
	// for efficiency studies
	validation_hists_["matchedRec_SimPt_RD"]->Fill(simpt[simtrack]);
	validation_hists_["matchedRec_SimPhi_RD"]->Fill(simphi[simtrack]);
	validation_hists_["matchedRec_SimEta_RD"]->Fill(simeta[simtrack]);
        
	// for fake rate studies
	validation_hists_["matchedRec_RecPt_RD"]->Fill(candpt);
	validation_hists_["matchedRec_RecPhi_RD"]->Fill(candphi);
	validation_hists_["matchedRec_RecEta_RD"]->Fill(candeta);
	
	// count the matched simtrack!
	associated_indices_found_RD[simtrack]++;
      } // end if block for simTrack with simtrack found
      else{ // reco track currently already matched sim track, don't count towards efficiency, but include reco info to not count it as fake
	// for fake rate studies
	validation_hists_["matchedRec_RecPt_RD"]->Fill(candpt);
	validation_hists_["matchedRec_RecPhi_RD"]->Fill(candphi);
	validation_hists_["matchedRec_RecEta_RD"]->Fill(candeta);
	// count the matched simtrack!
	associated_indices_found_RD[simtrack]++;
      } // end count duplicates
    } // end criterion if block
    //    }
    // Check to see if reco track has enough hits matched -- SD

    unsigned int denom_nHits_SD = 0;
    denom_nHits_SD = evt_sim_tracks[simtrack].nHits();    
    //    if (denom_nHits_RD >= 10){
    if (4*nHitsMatched >= 3*denom_nHits_SD){ // if association criterion is passed, save the info
      if (associated_indices_found_SD[simtrack] == 0){ // currently unmatched simtrack, count it towards efficiency 
	// for efficiency studies
	validation_hists_["matchedRec_SimPt_SD"]->Fill(simpt[simtrack]);
	validation_hists_["matchedRec_SimPhi_SD"]->Fill(simphi[simtrack]);
	validation_hists_["matchedRec_SimEta_SD"]->Fill(simeta[simtrack]);
	
	// for fake rate studies
	validation_hists_["matchedRec_RecPt_SD"]->Fill(candpt);
	validation_hists_["matchedRec_RecPhi_SD"]->Fill(candphi);
	validation_hists_["matchedRec_RecEta_SD"]->Fill(candeta);
	
	// count the matched simtrack!
	associated_indices_found_SD[simtrack]++;
      }
      else{ // currently already matched sim simtrack, don't count towards efficiency, but include reco info to not count it as fake
	// for fake rate studies
	validation_hists_["matchedRec_RecPt_SD"]->Fill(candpt);
	validation_hists_["matchedRec_RecPhi_SD"]->Fill(candphi);
	validation_hists_["matchedRec_RecEta_SD"]->Fill(candeta);
	// count the matched simtrack!
	associated_indices_found_SD[simtrack]++;
      }
    }
    //    }
  } // end loop over reco track candidate collection

  for (auto index = 0U; index < evt_sim_tracks.size(); ++index){ // loop over keys in map for duplicate tracks
    if (associated_indices_found_RD[index] > 1){ // check if map of RD associated indices has a duplicate
      // fill n-1 times e.g. if assoc_found == 3, means three tracks matched to one simTrack, i.e. TWO duplicates of one simTrack
      for (auto iduplicates = 0U; iduplicates < associated_indices_found_RD[index] - 1; ++iduplicates){
        validation_hists_["duplicateRec_SimPt_RD"]->Fill(simpt[index]);
        validation_hists_["duplicateRec_SimPhi_RD"]->Fill(simphi[index]);
        validation_hists_["duplicateRec_SimEta_RD"]->Fill(simeta[index]);
      }
    }
    if (associated_indices_found_SD[index] > 1){ // check if map of SD associated indices has a duplicate
      for (auto iduplicates = 0U; iduplicates < associated_indices_found_SD[index] - 1; ++iduplicates){
        validation_hists_["duplicateRec_SimPt_SD"]->Fill(simpt[index]);
        validation_hists_["duplicateRec_SimPhi_SD"]->Fill(simphi[index]);
	validation_hists_["duplicateRec_SimEta_SD"]->Fill(simeta[index]);
      }
    }
  } // end loop over index check for duplicates

  // for unassociated sim tracks
  for (auto index = 0U; index < evt_sim_tracks.size(); ++index){ // loop over keys in map for duplicate tracks
    if (associated_indices_found_RD[index] == 0){ // check if map of RD associated indices has a sim track without any rec tracks
      validation_hists_["missedRec_SimPt_RD"]->Fill(simpt[index]);
      validation_hists_["missedRec_SimPhi_RD"]->Fill(simphi[index]);
      validation_hists_["missedRec_SimEta_RD"]->Fill(simeta[index]);
      for (auto&& simhit : evt_sim_tracks[index].hitsVector()){
	validation_hists2_["missedRec_SimXY_RD"]->Fill(simhit.position()[0],simhit.position()[1]);
	validation_hists2_["missedRec_SimZR_RD"]->Fill(simhit.position()[2],simhit.r());
      }
    }
    if (associated_indices_found_SD[index] == 0){ // check if map of SD associated indices has a sim track without any rec track
      validation_hists_["missedRec_SimPt_SD"]->Fill(simpt[index]);
      validation_hists_["missedRec_SimPhi_SD"]->Fill(simphi[index]);
      validation_hists_["missedRec_SimEta_SD"]->Fill(simeta[index]);
      for (auto&& simhit : evt_sim_tracks[index].hitsVector()){
	validation_hists2_["missedRec_SimXY_SD"]->Fill(simhit.position()[0],simhit.position()[1]);
	validation_hists2_["missedRec_SimZR_SD"]->Fill(simhit.position()[2],simhit.r());
      }
    }
  } // end loop over index check for duplicates
}

void RootValidation::fillBuildHists(unsigned int layer, unsigned int branches, unsigned int cands)
{
  if (savetree_) {
    std::lock_guard<std::mutex> locker(glock_);
    layer_ = layer;
    branches_ = branches;
    cands_ = cands;
    tree_br_->Fill();
  }
}

void RootValidation::fillFitStateHists(const TrackState& simStateHit0, const TrackState& cfitStateHit0)
{
  // may get inconsistencies due to changes made before the fill is done.
  if (savetree_) {
    std::lock_guard<std::mutex> locker(glock_);
    simHit0_x=simStateHit0.parameters[0];
    simHit0_y=simStateHit0.parameters[1];
    simHit0_z=simStateHit0.parameters[2];
    simHit0_px=simStateHit0.parameters[3];
    simHit0_py=simStateHit0.parameters[4];
    simHit0_pz=simStateHit0.parameters[5];
    cfitHit0_x=cfitStateHit0.parameters[0];
    cfitHit0_y=cfitStateHit0.parameters[1];
    cfitHit0_z=cfitStateHit0.parameters[2];
    cfitHit0_px=cfitStateHit0.parameters[3];
    cfitHit0_py=cfitStateHit0.parameters[4];
    cfitHit0_pz=cfitStateHit0.parameters[5];
    cfitHit0_xe=sqrt(cfitStateHit0.errors[0][0]);
    cfitHit0_ye=sqrt(cfitStateHit0.errors[1][1]);
    cfitHit0_ze=sqrt(cfitStateHit0.errors[2][2]);
    cfitHit0_pxe=sqrt(cfitStateHit0.errors[3][3]);
    cfitHit0_pye=sqrt(cfitStateHit0.errors[4][4]);
    cfitHit0_pze=sqrt(cfitStateHit0.errors[5][5]);
  }
}

void RootValidation::fillFitHitHists(unsigned int hitid, const HitVec& mcInitHitVec, const MeasurementState& measState, const TrackState& propState, const TrackState& updatedState)
{
  if (savetree_){
    std::lock_guard<std::mutex> locker(glock_);
    MeasurementState initMeasState;
    for (auto&& mchit : mcInitHitVec){
      if(mchit.hitID() == hitid){
        initMeasState = mchit.measurementState();
        break;
      }
    }

    x_init   = initMeasState.parameters[0];
    x_mc     = measState.parameters[0];
    x_mcerr  = measState.errors[0][0]; // sigma^2 of x_mc (same with y,z)
    x_prop   = propState.parameters[0];
    x_perr   = propState.errors[0][0]; // sigma^2 of x_prop
    x_update = updatedState.parameters[0];
    x_uerr   = updatedState.errors[0][0]; // sigma^2 of x_update

    y_init   = initMeasState.parameters[1];
    y_mc     = measState.parameters[1];
    y_mcerr  = measState.errors[1][1];
    y_prop   = propState.parameters[1];
    y_perr   = propState.errors[1][1];
    y_update = updatedState.parameters[1];
    y_uerr   = updatedState.errors[1][1];

    z_init   = initMeasState.parameters[2];
    z_mc     = measState.parameters[2];
    z_mcerr  = measState.errors[2][2];
    z_prop   = propState.parameters[2];
    z_perr   = propState.errors[2][2];
    z_update = updatedState.parameters[2];
    z_uerr   = updatedState.errors[2][2];

    xy_mcerr = measState.errors[0][1];

    r_init   = sqrt( initMeasState.parameters[0]*initMeasState.parameters[0] +
                     initMeasState.parameters[1]*initMeasState.parameters[1] );
    r_mc     = sqrt( measState.parameters[0]*measState.parameters[0] + 
                     measState.parameters[1]*measState.parameters[1] ); 
    r_prop   = sqrt( propState.parameters[0]*propState.parameters[0] + 
                     propState.parameters[1]*propState.parameters[1] );
    r_update = sqrt( updatedState.parameters[0]*updatedState.parameters[0] + 
                     updatedState.parameters[1]*updatedState.parameters[1] );

    phi_init   = getPhi(initMeasState.parameters[0],initMeasState.parameters[1]);
    phi_mc     = getPhi(measState.parameters[0],measState.parameters[1]);
    phi_mcerr  = ( measState.errors[0][0]*measState.parameters[0]*measState.parameters[0] +
                   measState.errors[1][1]*measState.parameters[1]*measState.parameters[1] - 
                   measState.errors[0][1]*measState.parameters[0]*measState.parameters[1] - 
                   measState.errors[1][0]*measState.parameters[1]*measState.parameters[0] ) / (r_mc*r_mc); // sigma^2 of phi
    phi_prop   = getPhi(propState.parameters[0],propState.parameters[1]);
    phi_perr   = ( propState.errors[0][0]*propState.parameters[0]*propState.parameters[0] +
                   propState.errors[1][1]*propState.parameters[1]*propState.parameters[1] - 
                   propState.errors[0][1]*propState.parameters[0]*propState.parameters[1] - 
                   propState.errors[1][0]*propState.parameters[1]*propState.parameters[0] ) / (r_prop*r_prop); // sigma^2 of phi
    phi_update = getPhi(updatedState.parameters[0],updatedState.parameters[1]);
    phi_uerr   = ( updatedState.errors[0][0]*updatedState.parameters[0]*updatedState.parameters[0] +
                   updatedState.errors[1][1]*updatedState.parameters[1]*updatedState.parameters[1] - 
                   updatedState.errors[0][1]*updatedState.parameters[0]*updatedState.parameters[1] - 
                   updatedState.errors[1][0]*updatedState.parameters[1]*updatedState.parameters[0] ) / (r_update*r_update); // sigma^2 of phi 
    posTree_->Fill();
  }
}

void RootValidation::fillFitTrackHists(const TrackState& initState, const TrackState& updatedState)
{
  if (savetree_) {
    std::lock_guard<std::mutex> locker(glock_);
    pt_mc  = sqrt(initState.parameters[3]*initState.parameters[3]+initState.parameters[4]*initState.parameters[4]);
    pt_fit = sqrt(updatedState.parameters[3]*updatedState.parameters[3]+updatedState.parameters[4]*updatedState.parameters[4]);
    pt_err = sqrt( updatedState.errors[3][3]*updatedState.parameters[3]*updatedState.parameters[3] +
                   updatedState.errors[4][4]*updatedState.parameters[4]*updatedState.parameters[4] + 
                   2*updatedState.errors[3][4]*updatedState.parameters[3]*updatedState.parameters[4] )/pt_fit;
    fittree_->Fill();
   }
}

void RootValidation::saveHists() {
  if (savetree_) {
    std::lock_guard<std::mutex> locker(glock_);
    f_->cd();
    for(auto&& mapitr : validation_hists_){
      mapitr.second->Write();
    }
    for(auto&& mapitr : validation_hists2_){
      mapitr.second->Write();
    }
    f_->Write();
    f_->Close();
  }
}

void RootValidation::deleteHists() {
  std::lock_guard<std::mutex> locker(glock_);
  for(auto&& mapitr : validation_hists_) {
    delete (mapitr.second);
  }
  validation_hists_.clear();
  for(auto&& mapitr : validation_hists2_) {
    delete (mapitr.second);
  }
  validation_hists2_.clear();
}

void RootValidation::setupHists()
{
  validation_hists_["gen_trk_Pt"] = makeHist("h_gen_trk_Pt", "P_{T} of generated tracks", 30, 0, 15, "P_{T} [GeV]", "Gen Tracks");
  validation_hists_["gen_trk_Px"] = makeHist("h_gen_trk_Px", "P_{x} of generated tracks", 30, -15, 15, "P_{x} [GeV]", "Gen Tracks");
  validation_hists_["gen_trk_Py"] = makeHist("h_gen_trk_Py", "P_{y} of generated tracks", 30, -15, 15, "P_{y} [GeV]", "Gen Tracks");
  validation_hists_["gen_trk_Pz"] = makeHist("h_gen_trk_Pz", "P_{z} of generated tracks", 30, -20, 20, "P_{z} [GeV]", "Gen Tracks");
  validation_hists_["gen_trk_phi"] = makeHist("h_gen_trk_Phi", "phi of generated tracks from px/py", 20, -4, 4, "#phi", "Gen Tracks");
  validation_hists_["gen_trk_eta"] = makeHist("h_gen_trk_Eta", "eta of generated tracks", 100, -5, 5, "#eta", "Gen Tracks");
  validation_hists_["gen_trk_dPhi"] = makeHist("h_gen_trk_dPhi", "#Delta#phi between tracks", 20, 0, 4, "#Delta#phi", " Gen Tracks");
  validation_hists_["gen_trk_mindPhi"] = makeHist("h_gen_trk_mindPhi", "smallest #Delta#phi between tracks", 40, 0, 0.1, "#Delta#phi", "Gen Tracks");
  validation_hists_["gen_trk_dR"] = makeHist("h_gen_trk_dR", "#DeltaR between tracks", 20, 0, 4, "#Delta R", "Gen Tracks");
  validation_hists_["gen_trk_mindR"] = makeHist("h_gen_trk_mindR", "smallest #DeltaR between tracks", 40, 0, 0.5, "#Delta R", "Gen Tracks");
  validation_hists_["gen_hits_rad"] = makeHist("h_gen_hits_rad", "Radius of Hits",400,0,40,"Radius","Hits");
  validation_hists_["gen_hits_rad_lay3"] = makeHist("h_gen_hits_rad_lay3", "Radius of Hits in Layer 3",100,11.9,12.1,"Radius","Gen Hits");
  validation_hists_["gen_hits_cov00"] = makeHist("h_gen_hits_cov00", "Cov(X,X) for All Hits",1000,0.0000001,0.0001,"Covariance (cm^{2}","Gen Hits");
  validation_hists_["gen_hits_cov11"] = makeHist("h_gen_hits_cov11", "Cov(Y,Y) for All Hits",1000,0.0000001,0.0001,"Covariance (cm^{2}","Gen Hits");
  validation_hists_["gen_trk_nHits"]  = makeHist("h_gen_trk_nHits", "number of hits in simulated track", 11, -0.5,10.5, "# Hits per Sim Track", "Sim Tracks");

  // reco distributions

  validation_hists_["rec_trk_nHits"] = makeHist("h_rec_trk_nHits", "number of hits in reco track", 11, -0.5,10.5, "# Hits per Track Candidate", "Reco Tracks");
  validation_hists_["rec_trk_Pt"]    = makeHist("h_rec_trk_Pt", "P_{T} of reconstructed tracks", 30, 0, 15, "P_{T} [GeV]", "Reco Tracks");
  validation_hists_["rec_trk_phi"]   = makeHist("h_rec_trk_Phi", "phi of reconstructed tracks from px/py", 20, -4, 4, "#phi", "Reco Tracks");
  validation_hists_["rec_trk_eta"]   = makeHist("h_rec_trk_Eta", "eta of reconstructed tracks", 100, -5, 5, "#eta", "Reco Tracks");
  validation_hists_["rec_trk_dphi"] = makeHist("h_rec_trk_dphi", "dphi of rec tracks from y/x", 200, -0.2, 0.2, "#phi", "Reco Tracks");
  validation_hists_["rec_trk_chi2"] = makeHist("h_rec_trk_chi2", "chi2 of rec tracks", 100, 0, 100, "#chi^{2}", "Reco Tracks");

  // association plots for eff and fake rate studies

  // changed eta from 100, -5, 5 from 100, -1, 1

  validation_hists_["matchedRec_SimPt_RD"]  = makeHist("h_matchedRec_SimPt_RD", "Sim P_{T} of associated reco tracks (RecDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["matchedRec_SimPhi_RD"] = makeHist("h_matchedRec_SimPhi_RD", "Sim phi of associated reco tracks from px/py (RecDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["matchedRec_SimEta_RD"] = makeHist("h_matchedRec_SimEta_RD", "Sim eta of associated reco tracks (RecDenom)", 100, -5, 5, "#eta", "Tracks");
  
  validation_hists_["matchedRec_RecPt_RD"]  = makeHist("h_matchedRec_RecPt_RD", "Rec P_{T} of associated reco tracks (RecDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["matchedRec_RecPhi_RD"] = makeHist("h_matchedRec_RecPhi_RD", "Rec phi of associated reco tracks from px/py (RecDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["matchedRec_RecEta_RD"] = makeHist("h_matchedRec_RecEta_RD", "Rec eta of associated reco tracks (RecDenom)", 100, -5, 5, "#eta", "Tracks");

  validation_hists_["matchedRec_SimPt_SD"]  = makeHist("h_matchedRec_SimPt_SD", "Sim P_{T} of associated reco tracks (SimDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["matchedRec_SimPhi_SD"] = makeHist("h_matchedRec_SimPhi_SD", "Sim phi of associated reco tracks from px/py (SimDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["matchedRec_SimEta_SD"] = makeHist("h_matchedRec_SimEta_SD", "Sim eta of associated reco tracks (SimDenom)", 100, -5, 5, "#eta", "Tracks");
  
  validation_hists_["matchedRec_RecPt_SD"]  = makeHist("h_matchedRec_RecPt_SD", "Rec P_{T} of associated reco tracks (SimDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["matchedRec_RecPhi_SD"] = makeHist("h_matchedRec_RecPhi_SD", "Rec phi of associated reco tracks from px/py (SimDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["matchedRec_RecEta_SD"] = makeHist("h_matchedRec_RecEta_SD", "Rec eta of associated reco tracks (SimDenom)", 100, -5, 5, "#eta", "Tracks");

  // Duplicate Tracks Plots
  validation_hists_["duplicateRec_SimPt_RD"]  = makeHist("h_duplicateRec_SimPt_RD", "Sim P_{T} of duplicate reco tracks (RecDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["duplicateRec_SimPhi_RD"] = makeHist("h_duplicateRec_SimPhi_RD", "Sim phi of duplicate reco tracks from px/py (RecDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["duplicateRec_SimEta_RD"] = makeHist("h_duplicateRec_SimEta_RD", "Sim eta of duplicate reco tracks (RecDenom)", 100, -5, 5, "#eta", "Tracks");

  validation_hists_["duplicateRec_SimPt_SD"]  = makeHist("h_duplicateRec_SimPt_SD", "Sim P_{T} of duplicate reco tracks (SimDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["duplicateRec_SimPhi_SD"] = makeHist("h_duplicateRec_SimPhi_SD", "Sim phi of duplicate reco tracks from px/py (SimDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["duplicateRec_SimEta_SD"] = makeHist("h_duplicateRec_SimEta_SD", "Sim eta of duplicate reco tracks (SimDenom)", 100, -5, 5, "#eta", "Tracks");

  // missed tracks

  validation_hists_["missedRec_SimPt_RD"] = makeHist("h_missedRec_SimPt_RD", "Sim P_{T} of unassociated sim tracks (RecDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["missedRec_SimPhi_RD"] = makeHist("h_missedRec_SimPhi_RD", "Sim phi of unassociated sim tracks from px/py (RecDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["missedRec_SimEta_RD"] = makeHist("h_missedRec_SimEta_RD", "Sim eta of unassociated  tracks (RecDenom)", 100, -5, 5, "#eta", "Tracks");
  validation_hists2_["missedRec_SimXY_RD"] = makeHist2("h_missedRec_SimXY_RD","Sim Hit XY Positions of unassociated sim tracks (RecDenom)",1000,-50,50,1000,-50,50,"x","y");
  validation_hists2_["missedRec_SimZR_RD"] = makeHist2("h_missedRec_SimZR_RD","Sim Hit Z vs Rad Positions of unassociated sim tracks (RecDenom)",5000,-250,250,1000,0,50,"z","r");

  validation_hists_["missedRec_SimPt_SD"] = makeHist("h_missedRec_SimPt_SD", "Sim P_{T} of unassociated sim tracks (SimDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["missedRec_SimPhi_SD"] = makeHist("h_missedRec_SimPhi_SD", "Sim phi of unassociated sim tracks from px/py (SimDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["missedRec_SimEta_SD"] = makeHist("h_missedRec_SimEta_SD", "Sim eta of unassociated  tracks (SimDenom)", 100, -5, 5, "#eta", "Tracks");
  validation_hists2_["missedRec_SimXY_SD"] = makeHist2("h_missedRec_SimXY_SD","Sim Hit XY Positions of unassociated sim tracks (SimDenom)",1000,-50,50,1000,-50,50,"x","y");
  validation_hists2_["missedRec_SimZR_SD"] = makeHist2("h_missedRec_SimZR_SD","Sim Hit Z vs Rad Positions of unassociated sim tracks (SimDenom)",5000,-250,250,1000,0,50,"z","r");

  // Hists for seeds 

  // Hit Pairs

  validation_hists_["matchedSeedPair_SimPt"]  = makeHist("h_matchedSeedPair_SimPt", "Sim P_{T} of associated seed hit pairs", 30, 0, 15, "P_{T} [GeV]", "Seed Pairs");
  validation_hists_["matchedSeedPair_SimPhi"] = makeHist("h_matchedSeedPair_SimPhi", "Sim phi of associated seed hit pairs", 20, -4, 4, "#phi", "Seed Pairs");
  validation_hists_["matchedSeedPair_SimEta"] = makeHist("h_matchedSeedPair_SimEta", "Sim eta of associated seed hit pairs", 100, -5, 5, "#eta", "Seed Pairs");

  validation_hists_["extraSeedPair_SimPt"]  = makeHist("h_extraSeedPair_SimPt", "Sim P_{T} of N extra seed hit pairs", 30, 0, 15, "P_{T} [GeV]", "Seed Pairs");
  validation_hists_["extraSeedPair_SimPhi"] = makeHist("h_extraSeedPair_SimPhi", "Sim phi of N extra seed hit pairs", 20, -4, 4, "#phi", "Seed Pairs");
  validation_hists_["extraSeedPair_SimEta"] = makeHist("h_extraSeedPair_SimEta", "Sim eta of N extra seed hit pairs", 100, -5, 5, "#eta", "Seed Pairs");

  // Hit Triplets

  validation_hists_["matchedSeedTriplet_SimPt"]  = makeHist("h_matchedSeedTriplet_SimPt", "Sim P_{T} of associated seed hit triplets", 30, 0, 15, "P_{T} [GeV]", "Seed Triplets");
  validation_hists_["matchedSeedTriplet_SimPhi"] = makeHist("h_matchedSeedTriplet_SimPhi", "Sim phi of associated seed hit triplets", 20, -4, 4, "#phi", "Seed Triplets");
  validation_hists_["matchedSeedTriplet_SimEta"] = makeHist("h_matchedSeedTriplet_SimEta", "Sim eta of associated seed hit triplets", 100, -5, 5, "#eta", "Seed Triplets");

  validation_hists_["extraSeedTriplet_SimPt"]  = makeHist("h_extraSeedTriplet_SimPt", "Sim P_{T} of N extra seed hit triplets", 30, 0, 15, "P_{T} [GeV]", "Seed Triplets");
  validation_hists_["extraSeedTriplet_SimPhi"] = makeHist("h_extraSeedTriplet_SimPhi", "Sim phi of N extra seed hit triplets", 20, -4, 4, "#phi", "Seed Triplets");
  validation_hists_["extraSeedTriplet_SimEta"] = makeHist("h_extraSeedTriplet_SimEta", "Sim eta of N extra seed hit triplets", 100, -5, 5, "#eta", "Seed Triplets");

  validation_hists_["missed_phi_SimPt"]  = makeHist("h_missed_phi_SimPt", "Sim P_{T} of N missed phi hit triplets", 30, 0, 15, "P_{T} [GeV]", "Seed Triplets");
  validation_hists_["missed_phi_SimPhi"] = makeHist("h_missed_phi_SimPhi", "Sim phi of N missed phi hit triplets", 20, -4, 4, "#phi", "Seed Triplets");
  validation_hists_["missed_phi_SimEta"] = makeHist("h_missed_phi_SimEta", "Sim eta of N missed phi hit triplets", 100, -5, 5, "#eta", "Seed Triplets");

  validation_hists_["missed_eta_SimPt"]  = makeHist("h_missed_eta_SimPt", "Sim P_{T} of N missed eta hit triplets", 30, 0, 15, "P_{T} [GeV]", "Seed Triplets");
  validation_hists_["missed_eta_SimPhi"] = makeHist("h_missed_eta_SimPhi", "Sim phi of N missed eta hit triplets", 20, -4, 4, "#phi", "Seed Triplets");
  validation_hists_["missed_eta_SimEta"] = makeHist("h_missed_eta_SimEta", "Sim eta of N missed eta hit triplets", 100, -5, 5, "#eta", "Seed Triplets");

  // Filtered and fitted seed tracks
  
  validation_hists_["matchedSeed_SimPt"]  = makeHist("h_matchedSeed_SimPt", "Sim P_{T} of associated seed tracks ", 30, 0, 15, "P_{T} [GeV]", "Seed Tracks");
  validation_hists_["matchedSeed_SimPhi"] = makeHist("h_matchedSeed_SimPhi", "Sim phi of associated seed tracks", 20, -4, 4, "#phi", "Seed Tracks");
  validation_hists_["matchedSeed_SimEta"] = makeHist("h_matchedSeed_SimEta", "Sim eta of associated seed tracks", 100, -5, 5, "#eta", "Seed Tracks");

  validation_hists_["extraSeed_SimPt"]  = makeHist("h_extraSeed_SimPt", "Sim P_{T} of N extra seed tracks", 30, 0, 15, "P_{T} [GeV]", "Seed Tracks");
  validation_hists_["extraSeed_SimPhi"] = makeHist("h_extraSeed_SimPhi", "Sim phi of N extra seed tracks", 20, -4, 4, "#phi", "Seed Tracks");
  validation_hists_["extraSeed_SimEta"] = makeHist("h_extraSeed_SimEta", "Sim eta of N extra seed tracks", 100, -5, 5, "#eta", "Seed Tracks");

  validation_hists_["deta"] = makeHist("h_deta_third", "Eta hit pair sim - sim third hit", 100, -0.1, 0.1, "d#eta", "Tracks");

  validation_hists2_["detectorzr"] = makeHist2("h_detectorzr","zrdetector",5000,-250,250,1000,0,50,"z","r");
  validation_hists2_["detectorxy"] = makeHist2("h_detectorxy","xydetector",1000,-50,50,1000,-50,50,"x","y");
  validation_hists2_["inc_detector"] = makeHist2("h_incdetector","incdetector",800,-40,40,800,-40,40,"x","y");
  validation_hists_["inc_detector_phi"] = makeHist("h_incphi", "Sim phi ", 2000, -4, 4, "#phi", "Tracks");

  validation_hists_["matchedSeedTriplet_chi2fit"] = makeHist("h_matchedchi2fit","matched triplet chi2 fit",100,0,100,"chi2","Seed Triplets");
  validation_hists_["extraSeedTriplet_chi2fit"] = makeHist("h_extrachi2fit","extra triplet chi2 fit",100,0,100,"chi2","Seed Triplets");

  validation_hists_["conf_pt_pull"] = makeHist("h_conf_pt_pull","Conformal fit p_{T} pull for truth seeds",100,-10,10,"p_{T} pull","Truth seeds");
  validation_hists_["conf_invpt_res"] = makeHist("h_conf_invpt_res","Conformal fit 1/p_{T} res for truth seeds",400,-1.0,1.0,"p_{T} res","Truth seeds");
  validation_hists_["conf_phi_res"] = makeHist("h_conf_phi_res","Conformal fit #phi res for truth seeds",400,-0.25,0.25,"#phi res","Truth seeds");
  validation_hists_["conf_theta_res"] = makeHist("h_conf_theta_res","Conformal fit #theta res for truth seeds",400,-0.25,0.25,"#theta res","Truth seeds");
  validation_hists_["seed_pt_pull"] = makeHist("h_confKalman_pt_pull","Conformal-Kalman fit p_{T} pull for truth seeds",100,-10,10,"p_{T} pull","Truth seeds");
  validation_hists_["seed_pt_res"] = makeHist("h_confKalman_pt_res","Conformal-Kalman fit p_{T} res for truth seeds",100,-5.0,5.0,"p_{T} res","Truth seeds");
  validation_hists_["seedMC_pt_pull"] = makeHist("h_KalmanMC_pt_pull","Kalman MC fit p_{T} pull for truth seeds",100,-10,10,"p_{T} pull","Truth seeds");
  validation_hists_["seedMC_pt_res"] = makeHist("h_KalmanMC_pt_res","Kalman MC fit p_{T} res for truth seeds",100,-0.5,0.5,"p_{T} res","Truth seeds");
}


TH1F* RootValidation::makeHist(const std::string& name, const std::string& title,
  const int nbins, const double min, const double max,
  const std::string& xlabel, const std::string& ylabel)
{
  TH1F* tmp = new TH1F(name.c_str(), title.c_str(), nbins, min, max);
  tmp->SetDirectory(NULL); //user is now responsible for deleting hists
  tmp->GetXaxis()->SetTitle(xlabel.c_str());
  tmp->GetYaxis()->SetTitle(ylabel.c_str());
  return tmp;
}

TH2F* RootValidation::makeHist2(const std::string& name, const std::string& title,
  const int nxbins, const double xmin, const double xmax,
  const int nybins, const double ymin, const double ymax,
  const std::string& xlabel, const std::string& ylabel)
{
  TH2F* tmp = new TH2F(name.c_str(), title.c_str(), nxbins, xmin, xmax, nybins, ymin, ymax);
  tmp->SetDirectory(NULL); //user is now responsible for deleting hists
  tmp->GetXaxis()->SetTitle(xlabel.c_str());
  tmp->GetYaxis()->SetTitle(ylabel.c_str());
  return tmp;
}
#endif
