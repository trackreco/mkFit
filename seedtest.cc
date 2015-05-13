#include "seedtest.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Hit.h"
#include "Event.h"

void buildSeedsByMC(const TrackVec& evt_sim_tracks, TrackVec& evt_seed_tracks){
  for (unsigned int itrack=0;itrack<evt_sim_tracks.size();++itrack) {
    const Track& trk = evt_sim_tracks[itrack];
    const HitVec& hits = trk.hitsVector();
    HitVec seedhits;
    float chi2 = 0;
    TrackState updatedState = trk.state();
    for (auto ilayer=0U;ilayer<Config::nlayers_per_seed;++ilayer) {//seeds have first three layers as seeds
      Hit seed_hit = hits[ilayer]; // do this for now to make initHits element number line up with HitId number
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = seed_hit.measurementState();
      updatedState = updateParameters(propState, measState);
      seedhits.push_back(seed_hit);
      chi2 += computeChi2(propState,measState);
    }
    Track seed(updatedState,seedhits,0);//fixme chi2
    evt_seed_tracks.push_back(seed);
  }
}

#ifdef YUP

void buildSeedsByTriplets(const HitVec& evt_lay_hits, TrackVec& evt_seed_tracks){
  //create seeds (from sim tracks for now)
  bool debug(false);

  std::vector<HitVec> seed_pairs;

  const float dZ = 1.0; // 3.25 - 3.3
  const float innerrad = 4.0; // average radius of inner radius

  const float alphaBeta = 0.0520195; // 0.0458378 --> for d0 = .0025 cm

  for (unsigned int ihit=0;ihit<layerHits_[1].size();++ihit) { // 1 = second layer
    float outerrad  = layerHits_[1][ihit].r();
    float outerphi  = layerHits_[1][ihit].phi();
    float outerhitz = layerHits_[1][ihit].z();

    

    for (index_iter = cand_hit_idx.begin(); index_iter != cand_hit_idx.end(); ++index_iter){
      HitVec seed_pair;
      seed_pair.push_back(layerHits_[0][*index_iter]);
      seed_pair.push_back(layerHits_[1][ihit]);
      seed_pairs.push_back(seed_pair);
    }
  }
  


  if (debug){
    std::cout << "Hit Pairs" << std::endl;
    for (unsigned int i = 0; i < seed_pairs.size(); i++){
      HitVec tempseed = seed_pairs[i];
      printf("ilay0: %1u ilay1: %1u  \n",
	     tempseed[0].mcTrackID(),
	     tempseed[1].mcTrackID()
	     );
    }
    std::cout << std::endl;
  }

  const float dPhi     = 0.0458712; 
  const float thirdRad = 12.0;
  std::vector<HitVec> seed_triplets;
  std::vector<HitVec> filtered_seeds;
  TrackVec conformalTracks;
  HitVec missed_phis;
  HitVec missed_etas;
  std::vector<float> deta;

  for (auto&& seed_pair : seed_pairs){

    float linePhi = getPhi(seed_pair[1].position()[0] - seed_pair[0].position()[0], seed_pair[1].position()[1] - seed_pair[0].position()[1]);
    float thirdPhiMinus = 0.0;
    float thirdPhiPlus  = 0.0;  

    if (seed_pair[0].phi() < seed_pair[1].phi() ){
      thirdPhiMinus = normalizedPhi(linePhi - dPhi);
      thirdPhiPlus  = normalizedPhi(linePhi);
    }
    else{
      thirdPhiMinus = normalizedPhi(linePhi);
      thirdPhiPlus  = normalizedPhi(linePhi + dPhi);
    }

    
    //    if (seed_pair[0].mcTrackID() == seed_pair[1].mcTrackID()){
    //float thirdPhi = simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2].phi();
      //      printf("phi1: %7.5f phi2: %7.5f phi3: %7.5f phi3M: %7.5f phi3P: %7.5f \n",
      //  printf("phi1: %7.5f phi2: %7.5f phi3: %2u phi3M: %2u phi3P: %2u \n",
      //     seed_pair[0].phi()+TMath::Pi(),
      //	     seed_pair[1].phi()+TMath::Pi(),
	     //  thirdPhi+TMath::Pi(),
	     // thirdPhiMinus+TMath::Pi(),
	     // thirdPhiPlus+TMath::Pi()
      //   getPhiPartition(thirdPhi),
	       ///     getPhiPartition(thirdPhiMinus),
      //     getPhiPartition(thirdPhiPlus)
      //     );
    //    }
     
    unsigned int phiBinMinus = getPhiPartition(thirdPhiMinus);
    unsigned int phiBinPlus  = getPhiPartition(thirdPhiPlus);

    // for dz displacements -- straight line window
    
    float thirdZline = 2*seed_pair[1].position()[2]-seed_pair[0].position()[2];
    /*
    if (seed_pair[0].mcTrackID() == seed_pair[1].mcTrackID()){
      float thirdZ = simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2].position()[2];
      float thirdr = simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2].r();
      //      printf("z1: %09.5f z2: %09.5f z3: %09.5f zline: %09.5f \n \n",
      printf("z1: %09.5f z2: %09.5f z3: %2u zline: %2u \n \n",
      //seed_pair[0].position()[2],
	 //    seed_pair[1].position()[2],
	   //  thirdZ,
	   //  thirdZline
	     seed_pair[0].eta(),
	     seed_pair[1].eta(),
	     getEtaPartition(getEta(thirdr,thirdZ),etaDet),
	     getEtaPartition(getEta(thirdRad,thirdZline),etaDet)
	     );
    }
    */


    // calculate deta for single tracks

    if (seed_pair[0].mcTrackID() == seed_pair[1].mcTrackID()){
      deta.push_back(getEta(thirdRad,thirdZline) - getEta(thirdRad,simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2].position()[2]));
    }

    //deta = 0.06 for almost max efficiency

#ifdef ETASEG
    unsigned int etaBinThird = getEtaPartition(getEta(thirdRad,thirdZline),etaDet);
    unsigned int etaBinMinus = getEtaPartition(normalizedEta(getEta(thirdRad,thirdZline)-0.06),etaDet);
    unsigned int etaBinPlus  = getEtaPartition(normalizedEta(getEta(thirdRad,thirdZline)+0.06),etaDet);
#else
    //    unsigned int etaBinThird = 0;
    unsigned int etabinMinus = 0;
    unsigned int etabinPlus  = 0;
#endif    

    if (seed_pair[0].mcTrackID() == seed_pair[1].mcTrackID()){
      unsigned int thirdPhiPart = getPhiPartition(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2].phi());
      unsigned int thirdEtaPart = getEtaPartition(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2].eta(),etaDet);
      if (phiBinMinus<=phiBinPlus){
	if ( (thirdPhiPart < phiBinMinus) || (thirdPhiPart > phiBinPlus) ){
	  missed_phis.push_back(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2]);
    	}
      }
      else{
	if ( (thirdPhiPart > phiBinMinus) || (thirdPhiPart < phiBinPlus) ){
	  missed_phis.push_back(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2]);
    	}
      }
      if ( (thirdEtaPart < etaBinMinus) || (thirdEtaPart > etaBinPlus) ){
	missed_etas.push_back(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2]);
      }
    }

    std::vector<unsigned int> cand_hit_idx;
    std::vector<unsigned int>::iterator index_iter; // iterator for vector


    for (unsigned int ieta = etaBinMinus; ieta <= etaBinPlus; ++ieta){
      BinInfo binInfoMinus = lay_eta_phi_hit_idx_[2][ieta][int(phiBinMinus)];
      BinInfo binInfoPlus  = lay_eta_phi_hit_idx_[2][ieta][int(phiBinPlus)];
      
      
      // Branch here from wrapping
      if (phiBinMinus<=phiBinPlus){
	unsigned int firstIndex = binInfoMinus.first;
	unsigned int maxIndex   = binInfoPlus.first+binInfoPlus.second;
	
	for (unsigned int ihit  = firstIndex; ihit < maxIndex; ++ihit){
	  cand_hit_idx.push_back(ihit);
	}
      } 
      else { // loop wrap around end of array for phiBinMinus > phiBinPlus
	unsigned int firstIndex = binInfoMinus.first;
	unsigned int etaBinSize = lay_eta_phi_hit_idx_[2][ieta][62].first+lay_eta_phi_hit_idx_[2][ieta][62].second;
	
	for (unsigned int ihit  = firstIndex; ihit < etaBinSize; ++ihit){
	  cand_hit_idx.push_back(ihit);
	}
	
	unsigned int etaBinStart= lay_eta_phi_hit_idx_[2][ieta][0].first;
	unsigned int maxIndex   = binInfoPlus.first+binInfoPlus.second;
	
	for (unsigned int ihit  = etaBinStart; ihit < maxIndex; ++ihit){
	  cand_hit_idx.push_back(ihit);
	}
      }
    }
    
    for (index_iter = cand_hit_idx.begin(); index_iter != cand_hit_idx.end(); ++index_iter){
      HitVec seed_triplet;
      seed_triplet.push_back(seed_pair[0]);
      seed_triplet.push_back(seed_pair[1]);
      seed_triplet.push_back(layerHits_[2][*index_iter]);
      seed_triplets.push_back(seed_triplet);
    }
  }

  if (debug) {
    std::cout << "Hit Triplets" << std::endl;
    for (unsigned int i = 0; i < seed_triplets.size(); i++){
      HitVec tempseed = seed_triplets[i];
      printf("ilay0: %1u ilay1: %1u ilay2: %1u \n",
	     tempseed[0].mcTrackID(),
	     tempseed[1].mcTrackID(),
	     tempseed[2].mcTrackID()
	     );
    }
  }

  // Seed cleaning --> do some chi2 fit for RZ line then throw away with high chi2
  // choose ind = r, dep = z... could do total chi2 but, z errors 10*r
  
  // A = y-int, B = slope

  // res on z is set to be hitposerrZ*hitposerrZ
  const float varZ    = 0.1*0.1;
  const float invsig2 = 3.*(1./varZ);

  std::vector<float> chi2triplet;

  for (auto&& seed_triplet : seed_triplets){
    // first do fit for three hits
    float sumx2sig2 = 0;
    float sumxsig2  = 0;
    float sumysig2  = 0;
    float sumxysig2 = 0;
    float chi2fit   = 0;
    for (auto&& seedhit : seed_triplet) { 
      sumx2sig2 += ((seedhit.r())*(seedhit.r()) / varZ);
      sumxsig2  += (seedhit.r() / varZ);
      sumysig2  += (seedhit.position()[2] / varZ);
      sumxysig2 += ((seedhit.r())*(seedhit.position()[2]) / varZ);
    }
    float norm   = 1./ ((invsig2*sumx2sig2) - (sumxsig2*sumxsig2));
    float aParam = norm * ((sumx2sig2 * sumysig2) - (sumxsig2*sumxysig2));
    float bParam = norm * ((invsig2 * sumxysig2)  - (sumxsig2*sumysig2));
      
    //now perform chi2 on fit!

    for (auto&& seedhit : seed_triplet) { 
      chi2fit += pow((seedhit.position()[2] - aParam - (bParam * seedhit.r())),2) / varZ;
    }
    
    //    std::cout << "slope: " << bParam << " y-int: " << aParam << std::endl;

    //float slope = (seed_triplet[2].position()[2] - seed_triplet[0].position()[2]) / (seed_triplet[2].r() - seed_triplet[0].r());
    // float yint  = seed_triplet[2].position()[2]  - (slope*seed_triplet[2].r());

    //    std::cout << "slope: " << slope << " y-int: " << yint << std::endl;

    //    std::cout << "chi2fit: " << chi2fit << std::endl;
   
    chi2triplet.push_back(chi2fit);
    if (chi2fit<9.0){
      filtered_seeds.push_back(seed_triplet);
    }
  }

  // now perform kalman fit on seeds
  //  for (auto&& seed_triplet : seed_triplets){
  for(auto&& seed_triplet : filtered_seeds){
    int charge = 0;
    if (seed_triplet[1].phi() > seed_triplet[2].phi()){charge = 1;}
    else {charge = -1;}

    TrackState cfitStateHit0;
    conformalFit(seed_triplet[0],seed_triplet[1],seed_triplet[2],charge,cfitStateHit0);

    TrackState updatedState = cfitStateHit0;

    Track conformalTrack(updatedState,seed_triplet,0.);
    conformalTracks.push_back(conformalTrack);

    for (auto ilayer=0U;ilayer<nlayers_per_seed;++ilayer) {
      Hit seed_hit = seed_triplet[ilayer];
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = seed_hit.measurementState();
      updatedState = updateParameters(propState, measState);
    }
    Track seed(updatedState,seed_triplet,0.);//fixme chi2
    seedTracks_.push_back(seed);
  }
}
#endif
