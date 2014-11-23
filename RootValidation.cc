#include "RootValidation.h"
#ifndef NO_ROOT
float getPt(float px, float py) { return sqrt(px*px + py*py); }
float getPhi(float px, float py) { return std::atan2(py, px); }
float getEta(float px, float py, float pz) {
  float theta = atan2( getPt(px,py), pz );
  return -1. * log( tan(theta/2.) );
}
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

void RootValidation::fillSimHists(TrackVec& evt_sim_tracks)
{
  for(unsigned int isim_track = 0; isim_track < evt_sim_tracks.size(); ++isim_track){
    // float gen_trk_Pt = sqrt( (evt_sim_tracks[isim_track].momentum()[0]) * (evt_sim_tracks[isim_track].momentum()[0]) +
    //                                               (evt_sim_tracks[isim_track].momentum()[1]) * (evt_sim_tracks[isim_track].momentum()[1]) );
    // float gen_trk_theta = atan2( gen_trk_Pt, evt_sim_tracks[isim_track].momentum()[2] );
    // float gen_trk_eta = -1. * log( tan(gen_trk_theta / 2.) );
    // validation_hists_["gen_trk_Pt"]->Fill( gen_trk_Pt );
    // validation_hists_["gen_trk_Px"]->Fill( evt_sim_tracks[isim_track].momentum()[0] );
    // validation_hists_["gen_trk_Py"]->Fill( evt_sim_tracks[isim_track].momentum()[1] ); 
    // validation_hists_["gen_trk_Pz"]->Fill( evt_sim_tracks[isim_track].momentum()[2] ); 
    // validation_hists_["gen_trk_phi"]->Fill( std::atan2(evt_sim_tracks[isim_track].momentum()[1], evt_sim_tracks[isim_track].momentum()[0]) ); //phi=arctan(y/x), atan2 returns -pi,pi
    // validation_hists_["gen_trk_eta"]->Fill( gen_trk_eta );
    
    validation_hists_["gen_trk_Pt"]->Fill( getPt(evt_sim_tracks[isim_track].momentum()[0], evt_sim_tracks[isim_track].momentum()[1]) );
    validation_hists_["gen_trk_Px"]->Fill( evt_sim_tracks[isim_track].momentum()[0] );
    validation_hists_["gen_trk_Py"]->Fill( evt_sim_tracks[isim_track].momentum()[1] ); 
    validation_hists_["gen_trk_Pz"]->Fill( evt_sim_tracks[isim_track].momentum()[2] ); 
    validation_hists_["gen_trk_phi"]->Fill( getPhi(evt_sim_tracks[isim_track].momentum()[0], evt_sim_tracks[isim_track].momentum()[1]) );
    validation_hists_["gen_trk_eta"]->Fill( getEta(evt_sim_tracks[isim_track].momentum()[0], evt_sim_tracks[isim_track].momentum()[1], evt_sim_tracks[isim_track].momentum()[2]) );
    
    HitVec& hits = evt_sim_tracks[isim_track].hitsVector();
    for (unsigned int ihit = 0; ihit < hits.size(); ihit++){
      float rad = sqrt(hits[ihit].position()[0]*hits[ihit].position()[0] + hits[ihit].position()[1]*hits[ihit].position()[1]);
      validation_hists_["gen_hits_rad"]->Fill( rad );
      // Fill histo for layer 3
      if ( (rad > 11.0) && (rad < 13.0) ) {
        validation_hists_["gen_hits_rad_lay3"]->Fill( rad );
      }
      validation_hists_["gen_hits_cov00"]->Fill( hits[ihit].error()[0][0] );
      validation_hists_["gen_hits_cov11"]->Fill( hits[ihit].error()[1][1] );
    }
    
    float mindR = 999999;
    float mindPhi = 999999;
    for( unsigned int jsim_track = 0; jsim_track < evt_sim_tracks.size(); ++jsim_track ){
      if(jsim_track != isim_track) {
        float phii=getPhi(evt_sim_tracks[isim_track].momentum()[0], evt_sim_tracks[isim_track].momentum()[1]);
        float etai=getEta(evt_sim_tracks[isim_track].momentum()[0], evt_sim_tracks[isim_track].momentum()[1], evt_sim_tracks[isim_track].momentum()[2]);
        float phij=getPhi(evt_sim_tracks[jsim_track].momentum()[0], evt_sim_tracks[jsim_track].momentum()[1]);
        float etaj=getEta(evt_sim_tracks[jsim_track].momentum()[0], evt_sim_tracks[jsim_track].momentum()[1], evt_sim_tracks[jsim_track].momentum()[2]);
        mindR=std::min(mindR, deltaR(phii, etai, phij, etaj));
        mindPhi=std::min(mindPhi, deltaPhi(phii, phij));
        if(jsim_track > isim_track){
          validation_hists_["gen_trk_dR"]->Fill( deltaR(phii, etai, phij, etaj) );
          validation_hists_["gen_trk_dPhi"]->Fill( deltaPhi(phii, phij) );
        }
      }
    } 
    validation_hists_["gen_trk_mindR"]->Fill( mindR );
    validation_hists_["gen_trk_mindPhi"]->Fill( mindPhi );
    validation_hists_["gen_trk_nHits"]->Fill( evt_sim_tracks[isim_track].nHits());
  }
}

void RootValidation::fillCandidateHists(TrackVec& evt_track_candidates)
{
  //dump candidates
  for (unsigned int itkcand=0; itkcand<evt_track_candidates.size(); ++itkcand) {
    Track tkcand = evt_track_candidates[itkcand];
    validation_hists_["rec_trk_nHits"]->Fill(tkcand.nHits());
    validation_hists_["rec_trk_chi2"]->Fill(tkcand.chi2());
    validation_hists_["rec_trk_phi"]->Fill( getPhi(tkcand.momentum()[0], tkcand.momentum()[1]) );
    validation_hists_["rec_trk_Pt"]->Fill( getPt(tkcand.momentum()[0], tkcand.momentum()[1]) );
    validation_hists_["rec_trk_eta"]->Fill( getEta(tkcand.momentum()[0], tkcand.momentum()[1], tkcand.momentum()[2]) );
    if (savetree_) {
      tk_nhits_ = tkcand.nHits();
      tk_chi2_ = tkcand.chi2();
      buildtree_->Fill();
    }
  }
}

void RootValidation::fillAssociationHists(TrackVec& evt_track_candidates, TrackVec& evt_sim_tracks){
  //setup for assocation 
  std::vector<unsigned int> simIndices; // vector has elements equal to to its index
  std::vector<unsigned int>::iterator index_iter;
  std::map <unsigned int,unsigned int> associated_indices_found_RD; // key is simIndex, value is n times simTrack has index matched to associated reco track
  std::map <unsigned int,unsigned int> associated_indices_found_SD; // key is simIndex, value is n times simTrack has index matched to associated reco track

  MCHitInfoVec mcHitInfoVec; // container for mcInfo for all MC hits
  MCHitInfoVec::iterator mcHitInfoVec_iter;

  // Initialize simIndices vector and fill up allMCHitInfo vector --> to be used for matching
  for (unsigned int itrack=0; itrack<evt_sim_tracks.size(); itrack++){
    simIndices.push_back(itrack);
    HitVec simHits = evt_sim_tracks[itrack].hitsVector();  
    for (unsigned int ihit=0; ihit<simHits.size(); ++ihit){ // push all simHits simInfo into one giant vector
      mcHitInfoVec.push_back(simHits[ihit].mcHitInfo());
    }
  }

  for (unsigned int itkcand=0;itkcand<evt_track_candidates.size();++itkcand) {
    Track tkcand     = evt_track_candidates[itkcand];
    HitVec& candHits = tkcand.hitsVector();

    // Do hit matching HERE

    std::map <unsigned int,unsigned int> hit_matches_found; // key is simIndex, value is n times reco track hits sim Index matches simTrack index
    for (unsigned int ihit=0;ihit<candHits.size();++ihit) {//loop over hits in candidates
      for (mcHitInfoVec_iter = mcHitInfoVec.begin(); mcHitInfoVec_iter != mcHitInfoVec.end(); ++mcHitInfoVec_iter){
	if (candHits[ihit].mcHitInfo() == *mcHitInfoVec_iter){ 
	  hit_matches_found[candHits[ihit].mcIndex()]++; // count the number of hits in a reco track that match a hit index in sim tracks
	}// END IF check on matched indices
      } // end loop over sim indices for check
    } // end loop over layers in reco track
    
    // Check to see if reco track has enough hits matched -- RD

    unsigned int denom_nHits_RD = tkcand.nHits();
    for (index_iter = simIndices.begin(); index_iter != simIndices.end(); ++index_iter){ // association = n rec hits match sim hits / n rec hits in track > 75%
      if ( (float)hit_matches_found[*index_iter] / (float)denom_nHits_RD >= .75){ // if association criterion is passed, save the info
	if (associated_indices_found_RD[*index_iter] == 0){ // currently unmatched index, count it towards efficiency 
	  // for efficiency studies
	  validation_hists_["matchedRec_SimPt_RD"]->Fill(getPt(evt_sim_tracks[*index_iter].momentum()[0],evt_sim_tracks[*index_iter].momentum()[1]));
	  validation_hists_["matchedRec_SimPhi_RD"]->Fill(getPhi(evt_sim_tracks[*index_iter].momentum()[0],evt_sim_tracks[*index_iter].momentum()[1]));
	  validation_hists_["matchedRec_SimEta_RD"]->Fill(getEta(evt_sim_tracks[*index_iter].momentum()[0],evt_sim_tracks[*index_iter].momentum()[1],evt_sim_tracks[*index_iter].momentum()[2]));
	
	  // for fake rate studies
	  validation_hists_["matchedRec_RecPt_RD"]->Fill(getPt(tkcand.momentum()[0],tkcand.momentum()[1]));
	  validation_hists_["matchedRec_RecPhi_RD"]->Fill(getPhi(tkcand.momentum()[0],tkcand.momentum()[1]));
	  validation_hists_["matchedRec_RecEta_RD"]->Fill(getEta(tkcand.momentum()[0],tkcand.momentum()[1],tkcand.momentum()[2]));

	  // count the matched index!
	  associated_indices_found_RD[*index_iter]++;
	} // end if block for simTrack with index found
	else{ // reco track currently already matched sim track, don't count towards efficiency, but include reco info to not count it as fake
	  // for fake rate studies
	  validation_hists_["matchedRec_RecPt_RD"]->Fill(getPt(tkcand.momentum()[0],tkcand.momentum()[1]));
	  validation_hists_["matchedRec_RecPhi_RD"]->Fill(getPhi(tkcand.momentum()[0],tkcand.momentum()[1]));
	  validation_hists_["matchedRec_RecEta_RD"]->Fill(getEta(tkcand.momentum()[0],tkcand.momentum()[1],tkcand.momentum()[2]));

	  // count the matched index!
	  associated_indices_found_RD[*index_iter]++;
	} // end cound duplicates
	// remove sim index for next search to avoid repeating sim indices for matching
	break; // after one found per track cand, no need to look for more 
      } // end criterion if blco
    } // end loop over search for dominant index in rec eff with rec denom

    // Check to see if reco track has enough hits matched -- SD

    unsigned int denom_nHits_SD = 0;
    for (index_iter = simIndices.begin(); index_iter != simIndices.end(); ++index_iter){ // association = n rec hits match sim hits / n rec hits in track > 75%
      denom_nHits_SD = evt_sim_tracks[*index_iter].nHits();    
      if ( (float)hit_matches_found[*index_iter] / (float)denom_nHits_SD >= .75){ // if association criterion is passed, save the info
	if (associated_indices_found_SD[*index_iter] == 0){ // currently unmatched index, count it towards efficiency 
	  // for efficiency studies
	  validation_hists_["matchedRec_SimPt_SD"]->Fill(getPt(evt_sim_tracks[*index_iter].momentum()[0],evt_sim_tracks[*index_iter].momentum()[1]));
	  validation_hists_["matchedRec_SimPhi_SD"]->Fill(getPhi(evt_sim_tracks[*index_iter].momentum()[0],evt_sim_tracks[*index_iter].momentum()[1]));
	  validation_hists_["matchedRec_SimEta_SD"]->Fill(getEta(evt_sim_tracks[*index_iter].momentum()[0],evt_sim_tracks[*index_iter].momentum()[1],evt_sim_tracks[*index_iter].momentum()[2]));
	
	  // for fake rate studies
	  validation_hists_["matchedRec_RecPt_SD"]->Fill(getPt(tkcand.momentum()[0],tkcand.momentum()[1]));
	  validation_hists_["matchedRec_RecPhi_SD"]->Fill(getPhi(tkcand.momentum()[0],tkcand.momentum()[1]));
	  validation_hists_["matchedRec_RecEta_SD"]->Fill(getEta(tkcand.momentum()[0],tkcand.momentum()[1],tkcand.momentum()[2]));

	  // count the matched index!
	  associated_indices_found_SD[*index_iter]++;
	}
	else{ // currently already matched sim index, don't count towards efficiency, but include reco info to not count it as fake
	  // for fake rate studies
	  validation_hists_["matchedRec_RecPt_SD"]->Fill(getPt(tkcand.momentum()[0],tkcand.momentum()[1]));
	  validation_hists_["matchedRec_RecPhi_SD"]->Fill(getPhi(tkcand.momentum()[0],tkcand.momentum()[1]));
	  validation_hists_["matchedRec_RecEta_SD"]->Fill(getEta(tkcand.momentum()[0],tkcand.momentum()[1],tkcand.momentum()[2]));

	  // count the matched index!
	  associated_indices_found_SD[*index_iter]++;
	}
	// remove sim index for next search to avoid repeating sim indices for matching
	break; // after one found per track cand, no need to look for more 
      }
    } // end loop over search for dominant index in rec eff with rec denom
  } // end loop over reco track candidate collection

  for (index_iter = simIndices.begin(); index_iter != simIndices.end(); ++index_iter){ // loop over keys in map for duplicate tracks
    if (associated_indices_found_RD[*index_iter] > 1){ // check if map of RD associated indices has a duplicate
      // fill n-1 times e.g. if assoc_found == 3, means three tracks matched to one simTrack, i.e. TWO duplicates of one simTrack
      for (unsigned int iduplicates = 0; iduplicates < associated_indices_found_RD[*index_iter] - 1; ++iduplicates){
	validation_hists_["duplicateRec_SimPt_RD"]->Fill(getPt(evt_sim_tracks[*index_iter].momentum()[0],evt_sim_tracks[*index_iter].momentum()[1]));
	validation_hists_["duplicateRec_SimPhi_RD"]->Fill(getPhi(evt_sim_tracks[*index_iter].momentum()[0],evt_sim_tracks[*index_iter].momentum()[1]));
	validation_hists_["duplicateRec_SimEta_RD"]->Fill(getEta(evt_sim_tracks[*index_iter].momentum()[0],evt_sim_tracks[*index_iter].momentum()[1],evt_sim_tracks[*index_iter].momentum()[2]));
      }
    }
    if (associated_indices_found_SD[*index_iter] > 1){ // check if map of SD associated indices has a duplicate
      for (unsigned int iduplicates = 0; iduplicates < associated_indices_found_SD[*index_iter] - 1; ++iduplicates){
	validation_hists_["duplicateRec_SimPt_SD"]->Fill(getPt(evt_sim_tracks[*index_iter].momentum()[0],evt_sim_tracks[*index_iter].momentum()[1]));
	validation_hists_["duplicateRec_SimPhi_SD"]->Fill(getPhi(evt_sim_tracks[*index_iter].momentum()[0],evt_sim_tracks[*index_iter].momentum()[1]));
	validation_hists_["duplicateRec_SimEta_SD"]->Fill(getEta(evt_sim_tracks[*index_iter].momentum()[0],evt_sim_tracks[*index_iter].momentum()[1],evt_sim_tracks[*index_iter].momentum()[2]));
      }
    }
  } // end loop over index check for duplicates
}
void RootValidation::fillBuildHists(unsigned int layer, unsigned int branches, unsigned int cands)
{
  if (savetree_) {
    layer_ = layer;
    branches_ = branches;
    cands_ = cands;
    tree_br_->Fill();
  }
}

void RootValidation::fillFitStateHists(TrackState& simStateHit0, TrackState& cfitStateHit0)
{
  if (savetree_) {
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

void RootValidation::fillFitHitHists(MeasurementState& initMeasState, MeasurementState& measState, TrackState& propState, TrackState& updatedState)
{
  if (savetree_){
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

    phi_init   = atan2(initMeasState.parameters[1],initMeasState.parameters[0]);
    phi_mc     = atan2(measState.parameters[1],measState.parameters[0]);
    phi_mcerr  = ( measState.errors[0][0]*measState.parameters[0]*measState.parameters[0] +
                   measState.errors[1][1]*measState.parameters[1]*measState.parameters[1] - 
                   measState.errors[0][1]*measState.parameters[0]*measState.parameters[1] - 
                   measState.errors[1][0]*measState.parameters[1]*measState.parameters[0] ) / (r_mc*r_mc); // sigma^2 of phi
    phi_prop   = atan2(propState.parameters[1],propState.parameters[0]);
    phi_perr   = ( propState.errors[0][0]*propState.parameters[0]*propState.parameters[0] +
                   propState.errors[1][1]*propState.parameters[1]*propState.parameters[1] - 
                   propState.errors[0][1]*propState.parameters[0]*propState.parameters[1] - 
                   propState.errors[1][0]*propState.parameters[1]*propState.parameters[0] ) / (r_prop*r_prop); // sigma^2 of phi
    phi_update = atan2(updatedState.parameters[1],updatedState.parameters[0]);
    phi_uerr   = ( updatedState.errors[0][0]*updatedState.parameters[0]*updatedState.parameters[0] +
                   updatedState.errors[1][1]*updatedState.parameters[1]*updatedState.parameters[1] - 
                   updatedState.errors[0][1]*updatedState.parameters[0]*updatedState.parameters[1] - 
                   updatedState.errors[1][0]*updatedState.parameters[1]*updatedState.parameters[0] ) / (r_update*r_update); // sigma^2 of phi 
    posTree_->Fill();
  }
}

void RootValidation::fillFitTrackHists(TrackState& initState, TrackState& updatedState)
{
  if (savetree_) {
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
    f_->cd();
    auto mapitrend = validation_hists_.end();
        
    for(auto mapitr = validation_hists_.begin(); mapitr != mapitrend; ++mapitr){
      (mapitr->second)->Write();
    }
    f_->Write();
    f_->Close();
  }
}

void RootValidation::deleteHists() {
  auto mapitrend = validation_hists_.end();
 
  for(auto mapitr = validation_hists_.begin(); mapitr != mapitrend; ++mapitr) {
    delete (mapitr->second);
  }
  validation_hists_.clear();
}

void RootValidation::setupHists()
{
  validation_hists_["gen_trk_Pt"] = makeHist("h_gen_trk_Pt", "P_{T} of generated tracks", 30, 0, 15, "P_{T} [GeV]", "Gen Tracks");
  validation_hists_["gen_trk_Px"] = makeHist("h_gen_trk_Px", "P_{x} of generated tracks", 30, -15, 15, "P_{x} [GeV]", "Gen Tracks");
  validation_hists_["gen_trk_Py"] = makeHist("h_gen_trk_Py", "P_{y} of generated tracks", 30, -15, 15, "P_{y} [GeV]", "Gen Tracks");
  validation_hists_["gen_trk_Pz"] = makeHist("h_gen_trk_Pz", "P_{z} of generated tracks", 30, -20, 20, "P_{z} [GeV]", "Gen Tracks");
  validation_hists_["gen_trk_phi"] = makeHist("h_gen_trk_Phi", "phi of generated tracks from px/py", 20, -4, 4, "#phi", "Gen Tracks");
  validation_hists_["gen_trk_eta"] = makeHist("h_gen_trk_Eta", "eta of generated tracks", 40, -2, 2, "#eta", "Gen Tracks");
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
  validation_hists_["rec_trk_eta"]   = makeHist("h_rec_trk_Eta", "eta of reconstructed tracks", 40, -2, 2, "#eta", "Reco Tracks");
  validation_hists_["rec_trk_dphi"] = makeHist("h_rec_trk_dphi", "dphi of rec tracks from y/x", 200, -0.2, 0.2, "#phi", "Reco Tracks");
  validation_hists_["rec_trk_chi2"] = makeHist("h_rec_trk_chi2", "chi2 of rec tracks", 100, 0, 100, "#chi^{2}", "Reco Tracks");

  // association plots for eff and fake rate studies

  validation_hists_["matchedRec_SimPt_RD"]  = makeHist("h_matchedRec_SimPt_RD", "Sim P_{T} of N matched reco tracks (RecDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["matchedRec_SimPhi_RD"] = makeHist("h_matchedRec_SimPhi_RD", "Sim phi of N matched reco tracks from px/py (RecDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["matchedRec_SimEta_RD"] = makeHist("h_matchedRec_SimEta_RD", "Sim eta of N matched reco tracks (RecDenom)", 40, -2, 2, "#eta", "Tracks");
  
  validation_hists_["matchedRec_RecPt_RD"]  = makeHist("h_matchedRec_RecPt_RD", "Rec P_{T} of N matched reco tracks (RecDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["matchedRec_RecPhi_RD"] = makeHist("h_matchedRec_RecPhi_RD", "Rec phi of N matched reco tracks from px/py (RecDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["matchedRec_RecEta_RD"] = makeHist("h_matchedRec_RecEta_RD", "Rec eta of N matched reco tracks (RecDenom)", 40, -2, 2, "#eta", "Tracks");

  validation_hists_["matchedRec_SimPt_SD"]  = makeHist("h_matchedRec_SimPt_SD", "Sim P_{T} of N matched reco tracks (SimDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["matchedRec_SimPhi_SD"] = makeHist("h_matchedRec_SimPhi_SD", "Sim phi of N matched reco tracks from px/py (SimDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["matchedRec_SimEta_SD"] = makeHist("h_matchedRec_SimEta_SD", "Sim eta of N matched reco tracks (SimDenom)", 40, -2, 2, "#eta", "Tracks");
  
  validation_hists_["matchedRec_RecPt_SD"]  = makeHist("h_matchedRec_RecPt_SD", "Rec P_{T} of N matched reco tracks (SimDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["matchedRec_RecPhi_SD"] = makeHist("h_matchedRec_RecPhi_SD", "Rec phi of N matched reco tracks from px/py (SimDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["matchedRec_RecEta_SD"] = makeHist("h_matchedRec_RecEta_SD", "Rec eta of N matched reco tracks (SimDenom)", 40, -2, 2, "#eta", "Tracks");

  // Duplicate Tracks Plots
  validation_hists_["duplicateRec_SimPt_RD"]  = makeHist("h_duplicateRec_SimPt_RD", "Sim P_{T} of N duplicate reco tracks (RecDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["duplicateRec_SimPhi_RD"] = makeHist("h_duplicateRec_SimPhi_RD", "Sim phi of N duplicate reco tracks from px/py (RecDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["duplicateRec_SimEta_RD"] = makeHist("h_duplicateRec_SimEta_RD", "Sim eta of N duplicate reco tracks (RecDenom)", 40, -2, 2, "#eta", "Tracks");

  validation_hists_["duplicateRec_SimPt_SD"]  = makeHist("h_duplicateRec_SimPt_SD", "Sim P_{T} of N duplicate reco tracks (RecDenom)", 30, 0, 15, "P_{T} [GeV]", "Tracks");
  validation_hists_["duplicateRec_SimPhi_SD"] = makeHist("h_duplicateRec_SimPhi_SD", "Sim phi of N duplicate reco tracks from px/py (RecDenom)", 20, -4, 4, "#phi", "Tracks");
  validation_hists_["duplicateRec_SimEta_SD"] = makeHist("h_duplicateRec_SimEta_SD", "Sim eta of N duplicate reco tracks (RecDenom)", 40, -2, 2, "#eta", "Tracks");


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
#endif
