#include "seedtest.h"
#include "Hit.h"
#include "Event.h"
#include "ConformalUtils.h"
#include "KalmanUtils.h"
#include "Propagation.h"

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

void buildSeedsByRoadTriplets(const std::vector<HitVec>& evt_lay_hits, const BinInfoMap& segmentMap, TrackVec& evt_seed_tracks){
  //create seeds (from sim tracks for now)
  bool debug(false);

  std::vector<HitVec> hitTriplets;   // first will be pairs, then triplets, then filtered chi2 triplets, then Conf fit, then KF fit
  buildHitPairs(evt_lay_hits,segmentMap[0],hitTriplets); // pass only first layer map ... no need to pass the whole thing!

#ifdef DEBUG  
  if (debug){
    std::cout << "Hit Pairs" << std::endl;
    for(auto&& hitPair : hitTriplets){
      printf("ilay0: %1u ilay1: %1u  \n",
	     hitPair[0].mcTrackID(),
	     hitPair[1].mcTrackID()
	     );
    }
    std::cout << std::endl;
  }
#endif
  
  buildHitTriplets(evt_lay_hits,segmentMap[2],hitTriplets);

#ifdef DEBUG
  if (debug) {
    std::cout << "Hit Triplets" << std::endl;
    for(auto&& hitTriplet : hitTriplets){
      printf("ilay0: %1u ilay1: %1u ilay2: %1u \n",
	     hitTriplet[0].mcTrackID(),
	     hitTriplet[1].mcTrackID(),
	     hitTriplet[2].mcTrackID()
	     );
    }
  }
#endif
  
  filterHitTripletsByRZChi2(hitTriplets); // filter based on RZ chi2 cut
  buildSeedsFromTriplets(hitTriplets,evt_seed_tracks);
}

void buildHitPairs(const std::vector<HitVec>& evt_lay_hits, const BinInfoLayerMap& segLayMap, std::vector<HitVec>& hit_pairs){
  const float dZ = 1.0; // 3.25 - 3.3
  const float innerrad = 4.0; // average radius of inner radius
  const float alphaBeta = 0.0520195; // 0.0458378 --> for d0 = .0025 cm --> analytically derived

  for (unsigned int ihit=0;ihit<evt_lay_hits[1].size();++ihit) { // 1 = second layer
    const float outerhitz = evt_lay_hits[1][ihit].z(); // remember, layer[0] is first layer! --> second layer = [1]
    const float outerphi  = evt_lay_hits[1][ihit].phi();

#ifdef ETASEG
    const auto etaBinMinus = getEtaPartition(normalizedEta(getEta(innerrad,(outerhitz-dZ)/2.)));
    const auto etaBinPlus  = getEtaPartition(normalizedEta(getEta(innerrad,(outerhitz+dZ)/2.)));
#else
    const auto etaBinMinus = 0U;
    const auto etaBinPlus  = 0U;
#endif
    const auto phiBinMinus = getPhiPartition(normalizedPhi(outerphi - alphaBeta));
    const auto phiBinPlus  = getPhiPartition(normalizedPhi(outerphi + alphaBeta));

    hitIndices cand_hit_idx = getCandHitIndices(etaBinMinus,etaBinPlus,phiBinMinus,phiBinPlus,segLayMap);
    hitIdxIter index_iter;  

    for (index_iter = cand_hit_idx.begin(); index_iter != cand_hit_idx.end(); ++index_iter){
      HitVec hit_pair;
      hit_pair.push_back(evt_lay_hits[0][*index_iter]);
      hit_pair.push_back(evt_lay_hits[1][ihit]);
      hit_pairs.push_back(hit_pair);
    }
  }
}  

void buildHitTriplets(const std::vector<HitVec>& evt_lay_hits, const BinInfoLayerMap& segLayMap, std::vector<HitVec>& hit_triplets){
  const float dEta     = 0.6; // for almost max efficiency --> empirically derived
  const float dPhi     = 0.0458712; // numerically+semianalytically derived
  const float thirdRad = 12.0; // average third radius

  std::vector<HitVec>::iterator hitTripsIter;
  
  for (hitTripsIter = hit_triplets.begin(); hitTripsIter != hit_triplets.end(); ++hitTripsIter){
    HitVec hit_pair = (*hitTripsIter);
#ifdef ETASEG
    const float thirdZline = 2*hit_pair[1].z()-hit_pair[0].z(); // for dz displacements -- straight line window
    const auto etaBinMinus = getEtaPartition(normalizedEta(getEta(thirdRad,thirdZline)-dEta));
    const auto etaBinPlus  = getEtaPartition(normalizedEta(getEta(thirdRad,thirdZline)+dEta));
#else
    const auto etaBinMinus = 0U;
    const auto etaBinPlus  = 0U;
#endif    
    const float linePhi = getPhi(hit_pair[1].position()[0] - hit_pair[0].position()[0], hit_pair[1].position()[1] - hit_pair[0].position()[1]);
    float thirdPhiMinus = 0.0;
    float thirdPhiPlus  = 0.0;  
    if (hit_pair[0].phi() < hit_pair[1].phi() ){
      thirdPhiMinus = normalizedPhi(linePhi - dPhi);
      thirdPhiPlus  = normalizedPhi(linePhi);
    }
    else{
      thirdPhiMinus = normalizedPhi(linePhi);
      thirdPhiPlus  = normalizedPhi(linePhi + dPhi);
    }
    const auto phiBinMinus = getPhiPartition(thirdPhiMinus);
    const auto phiBinPlus  = getPhiPartition(thirdPhiPlus);

    hitIndices cand_hit_idx = getCandHitIndices(etaBinMinus,etaBinPlus,phiBinMinus,phiBinPlus,segLayMap);
    hitIdxIter index_iter;  
    
    for (index_iter = cand_hit_idx.begin(); index_iter != cand_hit_idx.end(); ++index_iter){
      (*hitTripsIter).push_back(evt_lay_hits[2][*index_iter]); // copy in matched hit to seed pair .. do this for with copies 
      if (index_iter != cand_hit_idx.end()-1){ // make n copies of the first seed pair to keep adding hits
	hitTripsIter = hit_triplets.emplace(hitTripsIter+1,hit_pair); 
      }
    }
  }
}

void filterHitTripletsByRZChi2(std::vector<HitVec>& hit_triplets){
  // Seed cleaning --> do some chi2 fit for RZ line then throw away with high chi2
  // choose ind = r, dep = z... could do total chi2 but, z errors 10*r
  // A = y-int, B = slope

  // res on z is set to be hitposerrZ*hitposerrZ
  const float varZ    = 0.1*0.1; // from simulation ... will need to change if that changes!
  const float invsig2 = 3.*(1./varZ);

  std::vector<HitVec>::iterator hitTripsIter;

  for (hitTripsIter = hit_triplets.begin(); hitTripsIter != hit_triplets.end();){
    // first do fit for three hits
    float sumx2sig2 = 0;
    float sumxsig2  = 0;
    float sumysig2  = 0;
    float sumxysig2 = 0;
    for (auto&& seedhit : (*hitTripsIter)) { 
      sumx2sig2 += ((seedhit.r())*(seedhit.r()) / varZ);
      sumxsig2  += (seedhit.r() / varZ);
      sumysig2  += (seedhit.z() / varZ);
      sumxysig2 += ((seedhit.r())*(seedhit.z()) / varZ);
    }
    const float norm   = 1./ ((invsig2*sumx2sig2) - (sumxsig2*sumxsig2));
    const float aParam = norm * ((sumx2sig2 * sumysig2) - (sumxsig2*sumxysig2));
    const float bParam = norm * ((invsig2 * sumxysig2)  - (sumxsig2*sumysig2));
 
    float chi2fit = 0;     
    for (auto&& seedhit : (*hitTripsIter)) { //now perform chi2 on fit!
      chi2fit += pow((seedhit.z() - aParam - (bParam * seedhit.r())),2) / varZ;
    }

    if (chi2fit>9.0){ // cut empirically derived
      hitTripsIter = hit_triplets.erase(hitTripsIter);
    }
    else{
      ++hitTripsIter;
    }
  }
}

void buildSeedsFromTriplets(const std::vector<HitVec> & hit_triplets, TrackVec & evt_seed_tracks){
  // now perform kalman fit on seeds --> first need initial parameters --> get from Conformal fitter!
  const bool cf = false; // use errors derived for seeding
  for(auto&& hit_triplet : hit_triplets){
    int charge = 0;
    if (hit_triplet[1].phi() > hit_triplet[2].phi()){charge = 1;}
    else {charge = -1;}

    TrackState updatedState;
    conformalFit(hit_triplet[0],hit_triplet[1],hit_triplet[2],charge,updatedState,bool(false),cf); // first bool is backward fit == false

    for (auto ilayer=0U;ilayer<Config::nlayers_per_seed;++ilayer) {
      Hit seed_hit = hit_triplet[ilayer];
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = seed_hit.measurementState();
      updatedState = updateParameters(propState, measState);
    }
    Track seed(updatedState,hit_triplet,0.);//fixme chi2
    evt_seed_tracks.push_back(seed);
  }
}


