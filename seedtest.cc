#include "seedtest.h"
#include "Hit.h"
#include "Event.h"
#include "ConformalUtils.h"
#include "KalmanUtils.h"
#include "Propagation.h"
//#define DEBUG
#include "Debug.h"

inline int getCharge(const Hit & hit0, const Hit & hit1, const Hit & hit2){
  return ((hit2.y()-hit0.y())*(hit2.x()-hit1.x())>(hit2.y()-hit1.y())*(hit2.x()-hit0.x())?1:-1);
}

inline float predz(const float z0, const float r0, const float z2, const float r2, const float predr) {
  return (predr-r0)*(z2-z0) / (r2-r0) + z0;
}

void buildSeedsByMC(const TrackVec& evt_sim_tracks, TrackVec& evt_seed_tracks, TrackExtraVec& evt_seed_extras, Event& ev){
  bool debug(true);
  bool cfseed(false); // use CF fit for initial estimate

  for (int itrack=0;itrack<evt_sim_tracks.size();++itrack) {
    const Track& trk = evt_sim_tracks[itrack];
    int   seedhits[Config::nLayers];
    float sumchi2 = 0;

    TrackState updatedState;
    if (cfseed) {
      conformalFit(ev.layerHits_[0][trk.getHitIdx(0)],ev.layerHits_[1][trk.getHitIdx(1)],ev.layerHits_[2][trk.getHitIdx(2)],trk.charge(),updatedState,false); 
    }
    else {
      updatedState = trk.state();
    }

    dprint("processing sim track # " << itrack << " par=" << trk.parameters());
    TSLayerPairVec updatedStates; // validation for position pulls

    for (auto ilayer=0;ilayer<Config::nlayers_per_seed;++ilayer) {//seeds have first three layers as seeds
      auto hitidx = trk.getHitIdx(ilayer);
      const Hit& seed_hit = ev.layerHits_[ilayer][hitidx];
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
      if (Config::super_debug) { ev.validation_.collectPropTSLayerVecInfo(ilayer,propState); }
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
	std::cout << "Seeding failed to propagate to layer: " << ilayer << " for sim track: " << itrack << std::endl;
	break;
      }
#endif
      const auto& measState = seed_hit.measurementState();
      const float chi2 = computeChi2(propState,measState); 
      //      sumchi2 += chi2; //--> could use this to make the chi2
      if (Config::super_debug) { ev.validation_.collectChi2LayerVecInfo(ilayer,chi2); }
      updatedState = updateParameters(propState, measState);
      if (Config::super_debug) { ev.validation_.collectUpTSLayerVecInfo(ilayer,updatedState); }
      seedhits[ilayer] = hitidx;

      updatedStates.push_back(std::make_pair(ilayer,updatedState)); // validation
    }
    ev.validation_.collectSeedTkTSLayerPairVecMapInfo(itrack,updatedStates); // use to collect position pull info

    Track seed(updatedState,0.,itrack,Config::nlayers_per_seed,seedhits);//fixme chi2 (could use sumchi2)
    dprint("created seed track # " << itrack << " par=" << seed.parameters());
    evt_seed_tracks.push_back(seed);
    evt_seed_extras.emplace_back(itrack); 
  }
}

void buildSeedsByRZFirstRPhiSecond(TrackVec& evt_seed_tracks, TrackExtraVec& evt_seed_extras,
				   const std::vector<HitVec>& evt_lay_hits, const BinInfoMap& segmentMap, Event& ev){

  // first loop over first layer hits, then over third layer hits, then over second layer hits.
  // reject first/third layer pairs if r-z line not within 2sigma of luminous region
  // see if second layer hit is within small residual of predicted second layer z pos from r-z line of pair (two sigma of smear)
  // assume radii equally spaced ... greatly reduces problem
  // otherwise will have to make this more general...

  int rz_path = 0; // 0,1,2: bigger number less brute force, however, trade off in more calculations

  TripletIdxVec hitTriplets;
  if (rz_path == 0) {
    for (int i = 0; i < evt_lay_hits[0].size(); i++){
      for (int j = 0; j < evt_lay_hits[2].size(); j++){
	if ( abs(((3.*evt_lay_hits[0][i].z()-evt_lay_hits[2][j].z())/2.)) > (Config::seed_z0cut)) {continue;}
	const float z1 = (evt_lay_hits[0][i].z() + evt_lay_hits[2][j].z()) / 2.;
	for (int k = 0; k < evt_lay_hits[1].size(); k++){
	  if (fabs(z1-evt_lay_hits[1][k].z()) < Config::lay2Zcut) {
	    TripletIdx triplet = {i,k,j};
	    hitTriplets.push_back(triplet);
	  }
	}
      }
    }
  }
  else if (rz_path == 1) { // don't cheat on r-z calculation ... but still explore all eta-phi bins
    for (int i = 0; i < evt_lay_hits[0].size(); i++){
      const float z0 = evt_lay_hits[0][i].z();
      const float r0 = evt_lay_hits[0][i].r();
      for (int j = 0; j < evt_lay_hits[2].size(); j++){
	const float z2 = evt_lay_hits[2][j].z();
	const float r2 = evt_lay_hits[2][j].r();
	if (fabs(predz(z0,r0,z2,r2,0.)) > (Config::seed_z0cut)) {continue;}
	const float z1 = predz(z0,r0,z2,r2,Config::fRadialSpacing*2.);
	for (int k = 0; k < evt_lay_hits[1].size(); k++) {
	  if (fabs(z1-evt_lay_hits[1][k].z()) < Config::lay2Zcut) { // five sigma inclusion
	    TripletIdx triplet = {i,k,j};
	    hitTriplets.push_back(triplet);
	  }
	}
      }
    }
  }
  else if (rz_path == 2){ // full implementation with eta-phi bins
    for (int i = 0; i < evt_lay_hits[0].size(); i++){
      const float z0 = evt_lay_hits[0][i].z();
      const float r0 = evt_lay_hits[0][i].r();
      const int etabin  = getEtaPartition(evt_lay_hits[0][i].eta());
      const int etabinM = (etabin-2)>0?etabin-2:((etabin-1)>0?etabin-1:etabin); 
      const int etabinP = (etabin+2)<Config::nEtaPart?etabin+2:((etabin+1)<Config::nEtaPart?etabin+1:etabin ); 
      std::vector<int> cand_lay2_indices = getCandHitIndices(etabinM,etabinP,0,Config::nPhiPart-1,segmentMap[1]);
      for (auto&& j : cand_lay2_indices){
	const float z2 = evt_lay_hits[2][j].z();
	const float r2 = evt_lay_hits[2][j].r();
	if (fabs(predz(z0,r0,z2,r2,0.)) > (Config::seed_z0cut)) {continue;}
	const float z1 = predz(z0,r0,z2,r2,Config::fRadialSpacing*2.);
	// since 95% of 2nd layer hits are within eta partition of pred, use just that one bin
	const int etabin  = getEtaPartition(getEta(Config::fRadialSpacing*2.,z1));
	const int etabinM = (etabin-1)>0?etabin-1:etabin; 
	const int etabinP = (etabin+1)<Config::nEtaPart?etabin+1:etabin; 
	std::vector<int> cand_lay1_indices = getCandHitIndices(etabinM,etabinP,0,Config::nPhiPart-1,segmentMap[1]);
	for (auto&& k : cand_lay1_indices){
	  if (fabs(z1-evt_lay_hits[1][k].z()) < Config::lay2Zcut) { // five sigma inclusion
	    TripletIdx triplet = {i,k,j};
	    hitTriplets.push_back(triplet);
	  }
	}
      }
    }
  }

  if (!Config::super_debug) {
    ev.validation_.fillSeedInfoTree(hitTriplets,ev);
  }

  TripletIdxVec filteredTriplets;
  filterHitTripletsByCircleParams(evt_lay_hits,hitTriplets,filteredTriplets);  
  //  filterHitTripletsByRZChi2(evt_lay_hits,hitTriplets,filteredTriplets); // filter based on RZ chi2 cut

  if (!Config::super_debug) {
    ev.validation_.fillSeedTree(hitTriplets,filteredTriplets,ev);
  }

  // turn triplets into track seeds by performing CF + KF fit  
  // buildSeedsFromTriplets(evt_lay_hits,filteredTriplets,evt_seed_tracks,evt_seed_extras,ev); 
}

void filterHitTripletsByCircleParams(const std::vector<HitVec>& evt_lay_hits, const TripletIdxVec& hit_triplets, TripletIdxVec& filtered_triplets){
  for (auto&& hit_triplet : hit_triplets){
    const float x0 = evt_lay_hits[0][hit_triplet[0]].x();
    const float y0 = evt_lay_hits[0][hit_triplet[0]].y();
    const float x1 = evt_lay_hits[1][hit_triplet[1]].x();
    const float y1 = evt_lay_hits[1][hit_triplet[1]].y();
    const float x2 = evt_lay_hits[2][hit_triplet[2]].x();
    const float y2 = evt_lay_hits[2][hit_triplet[2]].y();

    // now fit a circle, extract pT and d0 from center and radius
    const float mr = (y1-y0)/(x1-x0);
    const float mt = (y2-y1)/(x2-x1);
    const float a  = (mr*mt*(y2-y0) + mr*(x1+x2) - mt*(x0+x1))/(2.*(mr-mt));
    const float b  = -1.*(a-(x0+x1)/2.)/mr + (y0+y1)/2.;
    const float r  = getHypot(x0-a,y0-b);
    if ((r >= Config::maxCurvR) && (fabs(getHypot(a,b)-r) <= Config::seed_d0cut)) {
      filtered_triplets.push_back(hit_triplet);
    } // d0 cut 5mm, pT cut 0.5 GeV (radius of 0.5 GeV track)
  }
}

void buildSeedsByRoadTriplets(TrackVec& evt_seed_tracks, TrackExtraVec& evt_seed_extras, 
			      const std::vector<HitVec>& evt_lay_hits, const BinInfoMap& segmentMap, Event& ev){
  bool debug(false);
  bool curve=true; // choose between two ways to get from pairs to triplets -- true is closer to CMSSW
  // first will be pairs, then triplets, then filtered chi2 triplets, then Conf fit, then KF fit

  PairIdxVec hitPairs; 
  buildHitPairs(evt_lay_hits,segmentMap[0],hitPairs); // pass only first layer map ... no need to pass the whole thing!

#ifdef DEBUG  
  if (debug){
    dprint("Hit Pairs");
    for(auto&& hitPair : hitPairs){
      printf("ilay0: %1u ilay1: %1u  \n",
 	     ev.simHitsInfo_[evt_lay_hits[0][hitPair[0]].mcHitID()].mcTrackID(),
 	     ev.simHitsInfo_[evt_lay_hits[1][hitPair[1]].mcHitID()].mcTrackID()
	     );
    }
    dprint("");
  }

  PairIdxVec truthPairsFound;
  for(auto&& hitPair : hitPairs){
    if (ev.simHitsInfo_[evt_lay_hits[0][hitPair[0]].mcHitID()].mcTrackID() == ev.simHitsInfo_[evt_lay_hits[1][hitPair[1]].mcHitID()].mcTrackID()) {
      truthPairsFound.push_back(hitPair);
    }
  }
  std::swap(truthPairsFound,hitPairs);
#endif
  
  TripletIdxVec hitTriplets;
  if (curve) {
    buildHitTripletsCurve(evt_lay_hits,segmentMap[2],hitPairs,hitTriplets);
  }
  else {
    buildHitTripletsApproxWindow(evt_lay_hits,segmentMap[2],hitPairs,hitTriplets);
  }

#ifdef DEBUG
  if (debug){
    dprint("Hit Triplets");
    for(auto&& hitTriplet : hitTriplets){
      printf("ilay0: %1u ilay1: %1u ilay2: %1u \n",
 	     ev.simHitsInfo_[evt_lay_hits[0][hitTriplet[0]].mcHitID()].mcTrackID(),
 	     ev.simHitsInfo_[evt_lay_hits[1][hitTriplet[1]].mcHitID()].mcTrackID()
 	     ev.simHitsInfo_[evt_lay_hits[2][hitTriplet[2]].mcHitID()].mcTrackID()
	     );
    }
    dprint("");
  }
#endif

  TripletIdxVec filteredTriplets;
  TripletIdxVec filteredTriplets1;
  filterHitTripletsBySecondLayerZResidual(evt_lay_hits,hitTriplets,filteredTriplets1); // filter based on residual of 2nd hit and line in R-Z for 1st, 3rd hit
  filterHitTripletsByCircleParams(evt_lay_hits,filteredTriplets1,filteredTriplets); // filter based on circle fit to three hits (d0 and curvature)
  //  filterHitTripletsByRZChi2(evt_lay_hits,hitTriplets,filteredTriplets); // filter based on RZ chi2 cut
  if (!Config::super_debug) {
    ev.validation_.fillSeedTree(hitTriplets,filteredTriplets,ev);
  }
  buildSeedsFromTriplets(evt_lay_hits,filteredTriplets,evt_seed_tracks,evt_seed_extras,ev);
}

void buildHitPairs(const std::vector<HitVec>& evt_lay_hits, const BinInfoLayerMap& segLayMap, PairIdxVec& hit_pairs){
  // use only average radius of inner radius for calculation
  // lay12angdiff is a parameter for phi search window derived numerically from Mathematica... see one of my old talks

  for (unsigned int ihit=0;ihit<evt_lay_hits[1].size();++ihit) { // 1 = second layer
    const float outerhitz = evt_lay_hits[1][ihit].z(); // remember, layer[0] is first layer! --> second layer = [1]
    const float outerphi  = evt_lay_hits[1][ihit].phi();

#ifdef ETASEG   // z cut is the largest killer... up to 3 sigma is probably best
    const auto etaBinMinus = getEtaPartition(getEta(Config::fRadialSpacing,(outerhitz-Config::seed_z0cut)/2.));
    const auto etaBinPlus  = getEtaPartition(getEta(Config::fRadialSpacing,(outerhitz+Config::seed_z0cut)/2.));
#else
    const auto etaBinMinus = 0;
    const auto etaBinPlus  = 0;
#endif
    const auto phiBinMinus = getPhiPartition(normalizedPhi(outerphi - Config::lay12angdiff));
    const auto phiBinPlus  = getPhiPartition(normalizedPhi(outerphi + Config::lay12angdiff));

    std::vector<int> cand_hit_indices = getCandHitIndices(etaBinMinus,etaBinPlus,phiBinMinus,phiBinPlus,segLayMap);
    for (auto&& cand_hit_idx : cand_hit_indices){
      PairIdx hit_pair;
      hit_pair[0] = cand_hit_idx;
      hit_pair[1] = ihit;
      hit_pairs.push_back(hit_pair);
    }
  }
}  

void buildHitTripletsCurve(const std::vector<HitVec>& evt_lay_hits, const BinInfoLayerMap& segLayMap, 
			   const PairIdxVec& hit_pairs, TripletIdxVec& hit_triplets){
  const float lay3rad   = Config::fRadialSpacing * 3.0; // average third radius
  const float maxCurvR2 = Config::maxCurvR * Config::maxCurvR;

  for (auto&& hit_pair : hit_pairs){
    const Hit& hit0 = evt_lay_hits[0][hit_pair[0]];
    const Hit& hit1 = evt_lay_hits[1][hit_pair[1]];
    const float x0 = hit0.x(); const float y0 = hit0.y();
    const float x1 = hit1.x(); const float y1 = hit1.y();
    const float diffx2 = (x0-x1)*(x0-x1); 
    const float diffy2 = (y0-y1)*(y0-y1);

    // first solve for centers of max curvature circles fitting to the two points... (x[0,1]-a)^2 + (y[0,1]-b)^2 = r_max^2
    // then check if d0 of track is inside beamspot (maybe up to some sigma??) //fixme
    // if not, then force the track to be tangential to beamspot and through both points. --> system of 5 equations, yargh! //fixme
    // extrapolate x,y to third layer --> use intersection of two circles for each ... arbitrate solutions
    // rotate both sets of points, figure out which is the "less"er 
    // use getCandHitIndices (will handle phi-wrapping)

    const float quad = sqrtf((4*maxCurvR2 - diffx2 -diffy2) / (diffx2 + diffy2));

    // center of positive curved track
    const float apos = 0.5*((x0+x1)+(y0-y1)*quad);
    const float bpos = 0.5*((y0+y1)-(x0-x1)*quad);

    // center of negative curved track
    const float aneg = 0.5*((x0+x1)-(y0-y1)*quad);
    const float bneg = 0.5*((y0+y1)+(x0-x1)*quad);

    // now check if inside beamspot
    if (getHypot(apos,bpos)-Config::maxCurvR>getHypot(Config::beamspotX,Config::beamspotY)){
      // force pos inside beamspot...
    }
    if (getHypot(aneg,bneg)-Config::maxCurvR>getHypot(Config::beamspotX,Config::beamspotY)){
      // force neg inside beamspot...
    }

    // positive points of intersection with third layer
    float posx2 = 0., posy2 = 0;
    intersectThirdLayer(apos,bpos,x1,y1,posx2,posy2);
    const float posPhi = getPhi(posx2,posy2);

    // negative points of intersection with third layer
    float negx2 = 0., negy2 = 0;
    intersectThirdLayer(aneg,bneg,x1,y1,negx2,negy2);
    const float negPhi = getPhi(negx2,negy2);

#ifdef ETASEG
    const float thirdZline = 2*hit1.z()-hit0.z(); // for dz displacements -- straight line window
    const auto etaBinMinus = getEtaPartition(getEta(lay3rad,thirdZline)-Config::dEtaSeedTrip);
    const auto etaBinPlus  = getEtaPartition(getEta(lay3rad,thirdZline)+Config::dEtaSeedTrip);
#else
    const auto etaBinMinus = 0;
    const auto etaBinPlus  = 0;
#endif    
    const auto phiBinMinus = getPhiPartition(negPhi);
    const auto phiBinPlus  = getPhiPartition(posPhi);

#ifdef DEBUG
    const float lay2phi = evt_lay_hits[2][ev.simTracks_[ev.simHitsInfo_[hit0.mcHitID()].mcTrackID()].getHitIdx(2)].phi();
    dprint("lay0 phi: " << hit0.phi() << " lay1 phi: " << hit1.phi() << std::endl <<
	   "negPhi: " << negPhi << " lay2 phi: " << lay2phi << " posPhi: " << posPhi << std::endl <<
	   "binM: " << phiBinMinus << " phi2Bin: " << getPhiPartition(lay2phi) << " binP: " << phiBinPlus << std::endl);
#endif

    std::vector<int> cand_hit_indices = getCandHitIndices(etaBinMinus,etaBinPlus,phiBinMinus,phiBinPlus,segLayMap);
    for (auto&& cand_hit_idx : cand_hit_indices){
      TripletIdx hit_triplet;
      hit_triplet[0] = hit_pair[0];
      hit_triplet[1] = hit_pair[1];
      hit_triplet[2] = cand_hit_idx;      
      hit_triplets.push_back(hit_triplet);
    }
  }
}

void intersectThirdLayer(const float a, const float b, const float x1, const float y1, float& x2, float& y2){
  const float a2 = a*a; const float b2 = b*b; const float a2b2 = a2+b2;
  const float lay3rad2  = (Config::fRadialSpacing*Config::fRadialSpacing)*9.0; // average third radius squared
  const float maxCurvR2 = Config::maxCurvR * Config::maxCurvR;

  const float quad = sqrtf( 2*maxCurvR2*(a2b2+lay3rad2) - (a2b2-lay3rad2)*(a2b2-lay3rad2) - maxCurvR2*maxCurvR2 );
  const float pos[2] = { (a2*a + a*(b2+lay3rad2-maxCurvR2) - b*quad)/ a2b2 , (b2*b + b*(a2+lay3rad2-maxCurvR2) + a*quad)/ a2b2 };
  const float neg[2] = { (a2*a + a*(b2+lay3rad2-maxCurvR2) + b*quad)/ a2b2 , (b2*b + b*(a2+lay3rad2-maxCurvR2) - a*quad)/ a2b2 };

  // since we have two intersection points, arbitrate which one is closer to layer2 hit
  if (getHypot(pos[0]-x1,pos[1]-y1)<getHypot(neg[0]-x1,neg[1]-y1)) {
    x2 = pos[0];
    y2 = pos[1];
  }
  else {
    x2 = neg[0];
    y2 = neg[1];
  }
}

void buildHitTripletsApproxWindow(const std::vector<HitVec>& evt_lay_hits, const BinInfoLayerMap& segLayMap, 
				  const PairIdxVec& hit_pairs, TripletIdxVec& hit_triplets){
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
      TripletIdx hit_triplet;
      hit_triplet[0] = hit_pair[0];
      hit_triplet[1] = hit_pair[1];
      hit_triplet[2] = cand_hit_idx;      
      hit_triplets.push_back(hit_triplet);
    }
  }
}

void filterHitTripletsBySecondLayerZResidual(const std::vector<HitVec>& evt_lay_hits, const TripletIdxVec& hit_triplets, TripletIdxVec& filtered_triplets){
  for (auto&& hit_triplet : hit_triplets){
    //    const float z1 = predz(evt_lay_hits[0][hit_triplet[0]].z(),evt_lay_hits[0][hit_triplet[0]].r(),evt_lay_hits[2][hit_triplet[2]].z(),evt_lay_hits[2][hit_triplet[2]].r(),Config::fRadialSpacing*2.);
    const float z1 = (evt_lay_hits[0][hit_triplet[0]].z() + evt_lay_hits[2][hit_triplet[2]].z()) / 2.;
    if (fabs(z1-evt_lay_hits[1][hit_triplet[1]].z()) < Config::lay2Zcut) { // three sigma inclusion
      filtered_triplets.push_back(hit_triplet);
    }
  }
}
				     
void filterHitTripletsByRZChi2(const std::vector<HitVec>& evt_lay_hits, const TripletIdxVec& hit_triplets, 
			       TripletIdxVec& filtered_triplets){
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

void buildSeedsFromTriplets(const std::vector<HitVec>& evt_lay_hits, const TripletIdxVec& filtered_triplets, 
			    TrackVec& evt_seed_tracks, TrackExtraVec& evt_seed_extras, Event& ev){
  // now perform kalman fit on seeds --> first need initial parameters --> get from Conformal fitter!
  unsigned int seedID = 0;
  for(auto&& hit_triplet : filtered_triplets){
    int charge = getCharge(evt_lay_hits[0][hit_triplet[0]],evt_lay_hits[1][hit_triplet[1]],evt_lay_hits[2][hit_triplet[2]]);

    TrackState updatedState;
    conformalFit(evt_lay_hits[0][hit_triplet[0]],evt_lay_hits[1][hit_triplet[1]],evt_lay_hits[2][hit_triplet[2]],charge,updatedState,false); 

    // CF is bad at getting a good pT estimate, phi and theta are fine
    // "best case" config found in other studies is to set TS by hand:
    // x,y,z to (0,0,0)
    // exx,eyy,ezz to 4*(Config::beamspotXY/Z)^2 (2 sigma)
    // px,py,pz from CF the same
    // epxpx,epypy,epzpz set to 0.25*(px/py/pz)^2 (half sigma)

    ev.validation_.collectSeedTkCFMapInfo(seedID,updatedState);

    TSLayerPairVec updatedStates; // validation for position pulls
    float sumchi2 = 0;
    for (auto ilayer=0;ilayer<Config::nlayers_per_seed;++ilayer) {
      const Hit& seed_hit = evt_lay_hits[ilayer][hit_triplet[ilayer]];
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
      if (Config::super_debug) { ev.validation_.collectPropTSLayerVecInfo(ilayer,propState); }
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = seed_hit.measurementState();
      const float chi2 = computeChi2(propState,measState);
      sumchi2 += chi2;// --> could use this to make the chi2
      if (Config::super_debug) { ev.validation_.collectChi2LayerVecInfo(ilayer,chi2); } 
      updatedState = updateParameters(propState, measState);
      if (Config::super_debug) { ev.validation_.collectUpTSLayerVecInfo(ilayer,updatedState); }
      updatedStates.push_back(std::make_pair(ilayer,updatedState)); // validation
    }
    ev.validation_.collectSeedTkTSLayerPairVecMapInfo(seedID,updatedStates); // use to collect position pull info

    int hitIndices[3] = {hit_triplet[0],hit_triplet[1],hit_triplet[2]};
    Track seed(updatedState,sumchi2,seedID,Config::nlayers_per_seed,hitIndices);//fixme chi2
    evt_seed_tracks.push_back(seed);
    evt_seed_extras.emplace_back(seedID);
    seedID++; // increment dummy counter for seedID
  }
}

