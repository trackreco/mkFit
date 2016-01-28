#include "seedtest.h"
#include "Hit.h"
#include "Event.h"
#include "ConformalUtils.h"
#include "KalmanUtils.h"
#include "Propagation.h"
//#define DEBUG
#include "Debug.h"

void buildSeedsByMC(const TrackVec& evt_sim_tracks, TrackVec& evt_seed_tracks, TrackExtraVec& evt_seed_extras, Event& ev){
  bool debug(true);
  
  for (int itrack=0;itrack<evt_sim_tracks.size();++itrack) {
    const Track& trk = evt_sim_tracks[itrack];
    int   seedhits[Config::nLayers];
    float chi2 = 0;
    TrackState updatedState = trk.state();
    dprint("processing sim track # " << itrack << " par=" << trk.parameters());

    TSLayerPairVec updatedStates; // validation for position pulls

    for (auto ilayer=0;ilayer<Config::nlayers_per_seed;++ilayer) {//seeds have first three layers as seeds
      auto hitidx = trk.getHitIdx(ilayer);
      const Hit& seed_hit = ev.layerHits_[ilayer][hitidx];
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
	std::cout << "Seeding failed to propagate to layer: " << ilayer << " for sim track: " << itrack << std::endl;
	break;
      }
#endif
      const auto& measState = seed_hit.measurementState();
      chi2 += computeChi2(propState,measState); //--> could use this to make the chi2
      updatedState = updateParameters(propState, measState);
      seedhits[ilayer] = hitidx;

      updatedStates.push_back(std::make_pair(ilayer,updatedState)); // validation
    }
    ev.validation_.collectSeedTkTSLayerPairVecMapInfo(itrack,updatedStates); // use to collect position pull info

    Track seed(updatedState,0.,itrack,Config::nlayers_per_seed,seedhits);//fixme chi2
    dprint("created seed track # " << itrack << " par=" << seed.parameters());
    evt_seed_tracks.push_back(seed);
    evt_seed_extras.emplace_back(itrack); 
  }
}

void buildSeedsByRoadTriplets(TrackVec& evt_seed_tracks, TrackExtraVec& evt_seed_extras, 
			      const std::vector<HitVec>& evt_lay_hits, const BinInfoMap& segmentMap, Event& ev){
  bool debug(true);
  // first will be pairs, then triplets, then filtered chi2 triplets, then Conf fit, then KF fit

  std::vector<std::vector<int> > hitPairs; 
  buildHitPairs(evt_lay_hits,segmentMap[0],hitPairs); // pass only first layer map ... no need to pass the whole thing!

  

#ifdef DEBUG  
  if (debug){
    dprint("Hit Pairs");
    for(auto&& hitPair : hitPairs){
      printf("ilay0: %1u ilay1: %1u  \n",
 	     ev.simHitsInfo_[hitPair[0].mcHitID()].mcTrackID(),
 	     ev.simHitsInfo_[hitPair[1].mcHitID()].mcTrackID()
	     );
    }
    dprint("");
  }
#endif

  std::vector<std::vector<int> > hitTriplets;
  buildHitTriplets(evt_lay_hits,segmentMap[2],hitPairs,hitTriplets);

#ifdef DEBUG
  if (debug) {
       dprint("Hit Triplets");
       for(auto&& hitTriplet : hitTriplets){
	 printf("ilay0: %1u ilay1: %1u ilay2: %1u \n",
	     ev.simHitsInfo_[hitTriplet[0].mcHitID()].mcTrackID(),
	     ev.simHitsInfo_[hitTriplet[1].mcHitID()].mcTrackID()
	     ev.simHitsInfo_[hitTriplet[2].mcHitID()].mcTrackID()
	     );
    }
  }
#endif

  std::vector<std::vector<int> > filteredTriplets;
  filterHitTripletsByRZChi2(evt_lay_hits,hitTriplets,filteredTriplets); // filter based on RZ chi2 cut
  buildSeedsFromTriplets(evt_lay_hits,filteredTriplets,evt_seed_tracks,evt_seed_extras,ev);
}

void buildHitPairs(const std::vector<HitVec>& evt_lay_hits, const BinInfoLayerMap& segLayMap, std::vector<std::vector<int> >& hit_pairs){
  // use only average radius of inner radius for calculation
  // alphaBeta is a parameter for phi search window derived numerically from Mathematica... see one of my old talks

  for (unsigned int ihit=0;ihit<evt_lay_hits[1].size();++ihit) { // 1 = second layer
    const float outerhitz = evt_lay_hits[1][ihit].z(); // remember, layer[0] is first layer! --> second layer = [1]
    const float outerphi  = evt_lay_hits[1][ihit].phi();

#ifdef ETASEG   // z cut is the largest killer... up to 3 sigma is probably best
    const auto etaBinMinus = getEtaPartition(getEta(Config::fRadialSpacing,(outerhitz-Config::beamspotZ)/2.));
    const auto etaBinPlus  = getEtaPartition(getEta(Config::fRadialSpacing,(outerhitz+Config::beamspotZ)/2.));
#else
    const auto etaBinMinus = 0;
    const auto etaBinPlus  = 0;
#endif
    const auto phiBinMinus = getPhiPartition(normalizedPhi(outerphi - Config::alphaBeta));
    const auto phiBinPlus  = getPhiPartition(normalizedPhi(outerphi + Config::alphaBeta));

    std::vector<int> cand_hit_indices = getCandHitIndices(etaBinMinus,etaBinPlus,phiBinMinus,phiBinPlus,segLayMap);
    for (auto&& cand_hit_idx : cand_hit_indices){
      std::vector<int> hit_pair(2);
      hit_pair[0] = cand_hit_idx;
      hit_pair[1] = ihit;
      hit_pairs.push_back(hit_pair);
    }
  }
}  

void buildHitTriplets(const std::vector<HitVec>& evt_lay_hits, const BinInfoLayerMap& segLayMap, 
		      const std::vector<std::vector<int> >& hit_pairs, std::vector<std::vector<int> >& hit_triplets){
  const float thirdRad = Config::fRadialSpacing * 3.0; // average third radius

  for (auto&& hit_pair : hit_pairs){
    const Hit& hit0 = evt_lay_hits[0][hit_pair[0]];
    const Hit& hit1 = evt_lay_hits[1][hit_pair[1]];

#ifdef ETASEG
    const float thirdZline = 2*hit1.z()-hit0.z(); // for dz displacements -- straight line window
    const auto etaBinMinus = getEtaPartition(getEta(thirdRad,thirdZline)-Config::dEtaSeedTrip);
    const auto etaBinPlus  = getEtaPartition(getEta(thirdRad,thirdZline)+Config::dEtaSeedTrip);
#else
    const auto etaBinMinus = 0;
    const auto etaBinPlus  = 0;
#endif    
    const float linePhi = getPhi(hit1.position()[0] - hit0.position()[0], hit1.position()[1] - hit0.position()[1]);
    float thirdPhiMinus = 0.0;
    float thirdPhiPlus  = 0.0;  
    if (hit0.phi() < hit1.phi()){
      thirdPhiMinus = normalizedPhi(linePhi - Config::dPhiSeedTrip);
      thirdPhiPlus  = normalizedPhi(linePhi);
    }
    else{
      thirdPhiMinus = normalizedPhi(linePhi);
      thirdPhiPlus  = normalizedPhi(linePhi + Config::dPhiSeedTrip);
    }
    const auto phiBinMinus = getPhiPartition(thirdPhiMinus);
    const auto phiBinPlus  = getPhiPartition(thirdPhiPlus);

    std::vector<int> cand_hit_indices = getCandHitIndices(etaBinMinus,etaBinPlus,phiBinMinus,phiBinPlus,segLayMap);
    for (auto&& cand_hit_idx : cand_hit_indices){
      std::vector<int> hit_triplet(3);
      hit_triplet[0] = hit_pair[0];
      hit_triplet[1] = hit_pair[1];
      hit_triplet[2] = cand_hit_idx;      
      hit_triplets.push_back(hit_triplet);
    }
  }
}

void filterHitTripletsByRZChi2(const std::vector<HitVec>& evt_lay_hits, const std::vector<std::vector<int> >& hit_triplets, 
			       std::vector<std::vector<int> >& filtered_triplets){
  // Seed cleaning --> do some chi2 fit for RZ line then throw away with high chi2
  // choose ind = r, dep = z... could do total chi2 but, z errors 10*r
  // A = y-int, B = slope

  // res on z is set to be hitposerrZ*hitposerrZ
  const float invsig2 = 3.*(1./Config::varZ);

  for (auto&& hit_triplet : hit_triplets){
    // first do fit for three hits
    float sumx2sig2 = 0;
    float sumxsig2  = 0;
    float sumysig2  = 0;
    float sumxysig2 = 0;
    for (int i = 0; i < Config::nlayers_per_seed; i++) {
      const Hit& trip_hit = evt_lay_hits[i][hit_triplet[i]];

      sumx2sig2 += ((trip_hit.r())*(trip_hit.r()) / Config::varZ);
      sumxsig2  += (trip_hit.r() / Config::varZ);
      sumysig2  += (trip_hit.z() / Config::varZ);
      sumxysig2 += ((trip_hit.r())*(trip_hit.z()) / Config::varZ);
    }
    const float norm   = 1./ ((invsig2*sumx2sig2) - (sumxsig2*sumxsig2));
    const float aParam = norm * ((sumx2sig2 * sumysig2) - (sumxsig2*sumxysig2));
    const float bParam = norm * ((invsig2 * sumxysig2)  - (sumxsig2*sumysig2));
 
    float chi2fit = 0;     
    for (int i = 0; i < Config::nlayers_per_seed; i++) { //now perform chi2 on fit!
      const Hit& trip_hit = evt_lay_hits[i][hit_triplet[i]];
      chi2fit += pow((trip_hit.z() - aParam - (bParam * trip_hit.r())),2) / Config::varZ;
    }

    if (chi2fit<Config::chi2seedcut){ // cut empirically derived
      filtered_triplets.push_back(hit_triplet);
    }
  }
}

void buildSeedsFromTriplets(const std::vector<HitVec>& evt_lay_hits, const std::vector<std::vector<int> >& filtered_triplets, 
			    TrackVec& evt_seed_tracks, TrackExtraVec& evt_seed_extras, Event& ev){
  // now perform kalman fit on seeds --> first need initial parameters --> get from Conformal fitter!
  const bool fiterrs  = false; // use errors derived for seeding

  unsigned int seedID = 0;
  for(auto&& hit_triplet : filtered_triplets){
    int charge = getCharge(evt_lay_hits[0][hit_triplet[0]],evt_lay_hits[1][hit_triplet[1]],evt_lay_hits[2][hit_triplet[2]]);

    TrackState updatedState;
    conformalFit(evt_lay_hits[0][hit_triplet[0]],evt_lay_hits[1][hit_triplet[1]],evt_lay_hits[2][hit_triplet[2]],charge,updatedState,fiterrs); 

    // CF is bad at getting a good pT estimate, phi and theta are fine
    // "best case" config found in other studies is to set TS by hand:
    // x,y,z to (0,0,0)
    // exx,eyy,ezz to 4*(Config::beamspotXY/Z)^2 (2 sigma)
    // px,py,pz from CF the same
    // epxpx,epypy,epzpz set to 0.25*(px/py/pz)^2 (half sigma)

    ev.validation_.collectSeedTkCFMapInfo(seedID,updatedState);

    TSLayerPairVec updatedStates; // validation for position pulls
    float chi2 = 0;
    for (auto ilayer=0;ilayer<Config::nlayers_per_seed;++ilayer) {
      const Hit& seed_hit = evt_lay_hits[ilayer][hit_triplet[ilayer]];
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = seed_hit.measurementState();
      chi2 += computeChi2(propState,measState);// --> could use this to make the chi2
      updatedState = updateParameters(propState, measState);

      updatedStates.push_back(std::make_pair(ilayer,updatedState)); // validation
    }
    ev.validation_.collectSeedTkTSLayerPairVecMapInfo(seedID,updatedStates); // use to collect position pull info

    int hitIndices[3] = {hit_triplet[0],hit_triplet[1],hit_triplet[2]};
    Track seed(updatedState,chi2,seedID,Config::nlayers_per_seed,hitIndices);//fixme chi2
    evt_seed_tracks.push_back(seed);
    evt_seed_extras.emplace_back(seedID);
    seedID++; // increment dummy counter for seedID
  }
}

int getCharge(const Hit & hit0, const Hit & hit1, const Hit & hit2){
  float x0 = hit0.x(); float y0 = hit0.y();
  float x1 = hit1.x(); float y1 = hit1.y();
  float x2 = hit2.x(); float y2 = hit2.y();
  
  // first perform rotation of hits to first quadrant, placing phi(hit3) on pi/4, all others follow
  const float delphi = Config::PIOver4 - std::atan2(y2,x2);
  rotateHitPos(x0,y0,delphi);
  rotateHitPos(x1,y1,delphi);
  rotateHitPos(x2,y2,delphi);

  dprint("x0: " << x0 << " y0: " << y0);
  dprint("x1: " << x1 << " y1: " << y1);
  dprint("x2: " << x2 << " y2: " << y2);

  const float x      = extrapolateToInterceptX(x0,y0,x1,y1,x2,y2);
  const int   charge = (x1<x?1:-1);
  return charge;
}

void rotateHitPos(float& x, float& y, const float delphi){
  // perform rotation
  const float xprime = x*cos(delphi)-y*sin(delphi);
  const float yprime = x*sin(delphi)+y*cos(delphi);
  
  // set variables
  x = xprime;
  y = yprime;
}

float extrapolateToInterceptX(const float x0, const float y0, const float x1, const float y1, const float x2, const float y2){
  // line --> y = mx + b, m = slope, b = y-int
  const float m = (y2 - y0) / (x2 - x0);
  const float b = y2 - (m * x2);

  // need to solve for the intercept of a straight line with a circle
  // i.e. y = mx + b, and x^2 + y^2 = r^2, where we solve for x to get extrapolated point
  // (mx + b)^2 = r^2 - x^2 ...  
  // this amounts to solving the quadratic equation for x.
  // x = (-b +/- sqrt(b^2 - 4*a*c))/ 2a --> a = 1 + m^2, b = 2mb, c = b^2 - r^2

  const float r = sqrtf(getRad2(x1,y1)); //Config::fRadialSpacing * 2.0; // (2nd layer)

  const float det    = sqrtf(4*m*m*b*b - 4*(1+m*m)*(b*b - r*r));
  const float x_up   = (-2.0*m*b + det) / (2*(1+m*m));
  const float x_down = (-2.0*m*b - det) / (2*(1+m*m));
  const float y_up   = m*x_up   + b;
  const float y_down = m*x_down + b;

  // need to arbitrate two answers : select one which has a distance --> use outer layer distance
  const float dist_up   = sqrtf((x_up  -x2)*(x_up  -x2) + (y_up  -y2)*(y_up  -y2));
  const float dist_down = sqrtf((x_down-x2)*(x_down-x2) + (y_down-y2)*(y_down-y2));

  // since layer radius ~constant, one dist should be very close, while the other is off by x5
  return ( (dist_up < dist_down) ? x_up : x_down ); 
}
