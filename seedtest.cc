#include "seedtest.h"
#include "Hit.h"
#include "Event.h"
#include "ConformalUtils.h"
#include "KalmanUtils.h"
#include "Propagation.h"
//#define DEBUG
#include "Debug.h"

inline float predz(const float z0, const float r0, const float z2, const float r2, const float predr) {
  return (predr-r0)*(z2-z0) / (r2-r0) + z0;
}

void buildSeedsByMC(const TrackVec& evt_sim_tracks, TrackVec& evt_seed_tracks, TrackExtraVec& evt_seed_extras, Event& ev)
{
  fprintf(stderr, "__FILE__::__LINE__ Needs fixing for B/E support, search for XXMT4K\n");
  exit(1);

#ifdef DEBUG
  bool debug(true);
#endif

  const PropagationFlags pflags(PF_none);

  for (int itrack=0;itrack<evt_sim_tracks.size();++itrack)
  {
    const Track& trk = evt_sim_tracks[itrack];
    int   seedhits[Config::nLayers];
    //float sumchi2 = 0;

    TrackState updatedState;
    if (Config::cf_seeding) {
      conformalFit(ev.layerHits_[0][trk.getHitIdx(0)],ev.layerHits_[1][trk.getHitIdx(1)],ev.layerHits_[2][trk.getHitIdx(2)],updatedState,false); 
      updatedState.charge = trk.charge();
    }
    else {
      updatedState = trk.state();
    }

    dprint("processing sim track # " << itrack << " par=" << trk.parameters());
    TSLayerPairVec updatedStates; // validation for position pulls

    // use first layers as seed hits
    for (auto ilayer=0;ilayer<Config::nlayers_per_seed;++ilayer)
    {
      auto hitidx = trk.getHitIdx(ilayer);
      const Hit& seed_hit = ev.layerHits_[ilayer][hitidx];
      TrackState propState = propagateHelixToR(updatedState, seed_hit.r(), pflags);
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
	std::cout << "Seeding failed to propagate to layer: " << ilayer << " for sim track: " << itrack << std::endl;
	break;
      }
#endif
      const auto& measState = seed_hit.measurementState();
      //const float chi2 = computeChi2(propState,measState); 
      //      sumchi2 += chi2; //--> could use this to make the chi2
      updatedState = updateParameters(propState, measState);
      seedhits[ilayer] = hitidx;

      updatedStates.push_back(std::make_pair(ilayer,updatedState)); // validation
    }

    // XXMT4K: Here will need layer indices, too now.
    // XX Track seed(updatedState,0.0f,itrack,Config::nlayers_per_seed,seedhits);//fixme chi2 (could use sumchi2)
    // XX dprint("created seed track # " << itrack << " par=" << seed.parameters());
    // XX evt_seed_tracks.push_back(seed);
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
	if ( std::abs(((3.*evt_lay_hits[0][i].z()-evt_lay_hits[2][j].z())/2.)) > (Config::seed_z0cut)) {continue;}
	const float z1 = (evt_lay_hits[0][i].z() + evt_lay_hits[2][j].z()) / 2.;
	for (int k = 0; k < evt_lay_hits[1].size(); k++){
	  if (std::abs(z1-evt_lay_hits[1][k].z()) < Config::seed_z1cut) {
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
	if (std::abs(predz(z0,r0,z2,r2,0.0f)) > (Config::seed_z0cut)) {continue;}
	const float z1 = predz(z0,r0,z2,r2,Config::fRadialSpacing*2.);
	for (int k = 0; k < evt_lay_hits[1].size(); k++) {
	  if (std::abs(z1-evt_lay_hits[1][k].z()) < Config::seed_z1cut) { // five sigma inclusion
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
	if (std::abs(predz(z0,r0,z2,r2,0.0f)) > (Config::seed_z0cut)) {continue;}
	const float z1 = predz(z0,r0,z2,r2,Config::fRadialSpacing*2.);
	// since 95% of 2nd layer hits are within eta partition of pred, use just that one bin
	const int etabin  = getEtaPartition(getEta(Config::fRadialSpacing*2.,z1));
	const int etabinM = (etabin-1)>0?etabin-1:etabin; 
	const int etabinP = (etabin+1)<Config::nEtaPart?etabin+1:etabin; 
	std::vector<int> cand_lay1_indices = getCandHitIndices(etabinM,etabinP,0,Config::nPhiPart-1,segmentMap[1]);
	for (auto&& k : cand_lay1_indices){
	  if (std::abs(z1-evt_lay_hits[1][k].z()) < Config::seed_z1cut) { // five sigma inclusion
	    TripletIdx triplet = {i,k,j};
	    hitTriplets.push_back(triplet);
	  }
	}
      }
    }
  }

  TripletIdxVec filteredTriplets;
  filterHitTripletsByCircleParams(evt_lay_hits,hitTriplets,filteredTriplets);  
  //  filterHitTripletsByRZChi2(evt_lay_hits,hitTriplets,filteredTriplets); // filter based on RZ chi2 cut
  // turn triplets into track seeds by performing CF + KF fit  
  // buildSeedsFromTriplets(evt_lay_hits,filteredTriplets,evt_seed_tracks,evt_seed_extras,ev); 
}

void buildSeedsByRoadTriplets(TrackVec& evt_seed_tracks, TrackExtraVec& evt_seed_extras, 
			      const std::vector<HitVec>& evt_lay_hits, const BinInfoMap& segmentMap, Event& ev){
#ifdef DEBUG
  bool debug(false);
#endif
  bool curve=true; // choose between two ways to get from pairs to triplets -- true is closer to CMSSW
  // first will be pairs, then triplets, then filtered chi2 triplets, then Conf fit, then KF fit

  PairIdxVec hitPairs; 
  buildHitPairs(evt_lay_hits,segmentMap[0],hitPairs); // pass only first layer map ... no need to pass the whole thing!

#ifdef DEBUG  
  if (debug){
    dprint("Hit Pairs");
    for(auto&& hitPair : hitPairs){
      dprintf("ilay0: %1u ilay1: %1u  \n",
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
      dprintf("ilay0: %1u ilay1: %1u ilay2: %1u \n",
 	     ev.simHitsInfo_[evt_lay_hits[0][hitTriplet[0]].mcHitID()].mcTrackID(),
 	     ev.simHitsInfo_[evt_lay_hits[1][hitTriplet[1]].mcHitID()].mcTrackID(),
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
  buildSeedsFromTriplets(evt_lay_hits,filteredTriplets,evt_seed_tracks,evt_seed_extras,ev);
}

void buildSeedsByRoadSearch(TrackVec& evt_seed_tracks, TrackExtraVec& evt_seed_extras, 
			    const std::vector<HitVec>& evt_lay_hits, const BinInfoMap& segmentMap, Event& ev)
{
  fprintf(stderr, "__FILE__::__LINE__ Needs fixing for B/E support, search for XXMT4K\n");
  exit(1);

  // use this to initialize tracks
  const TrackState dummystate;

  const float lay2_r    = Config::fRadialSpacing * 3.0; // average third layer radius
  const float maxCurvR2 = Config::maxCurvR * Config::maxCurvR;

  int seedID = 0;
  for (int ihit1 = 0; ihit1 < evt_lay_hits[1].size(); ++ihit1) { // 0 = first layer, 1 = second layer, 2 = third layer
    const Hit & hit1     = evt_lay_hits[1][ihit1];
    const float hit1_z   = evt_lay_hits[1][ihit1].z();
    const float hit1_phi = evt_lay_hits[1][ihit1].phi();

    const auto lay0_etaBinM = getEtaPartition(getEta(Config::fRadialSpacing,(hit1_z-Config::seed_z0cut)/2.));
    const auto lay0_etaBinP = getEtaPartition(getEta(Config::fRadialSpacing,(hit1_z+Config::seed_z0cut)/2.));
    const auto lay0_phiBinM = getPhiPartition(normalizedPhi(hit1_phi - Config::lay01angdiff));
    const auto lay0_phiBinP = getPhiPartition(normalizedPhi(hit1_phi + Config::lay01angdiff));

    std::vector<int> cand_hit0_indices = getCandHitIndices(lay0_etaBinM,lay0_etaBinP,lay0_phiBinM,lay0_phiBinP,segmentMap[0]);
    for (auto&& ihit0 : cand_hit0_indices){

      const Hit & hit0     = evt_lay_hits[0][ihit0];
      const float hit0_z   = hit0.z();
      const float hit0_x   = hit0.x(); 
      const float hit0_y   = hit0.y();
      const float hit1_x   = hit1.x(); 
      const float hit1_y   = hit1.y();
      const float hit01_r2 = getRad2(hit0_x-hit1_x,hit0_y-hit1_y);

      const float quad = std::sqrt((4*maxCurvR2 - hit01_r2) / hit01_r2);
    
      // center of negative curved track
      const float aneg = 0.5*((hit0_x+hit1_x)-(hit0_y-hit1_y)*quad);
      const float bneg = 0.5*((hit0_y+hit1_y)+(hit0_x-hit1_x)*quad);

      // negative points of intersection with third layer
      float lay2_negx = 0.0f, lay2_negy = 0.0f;
      intersectThirdLayer(aneg,bneg,hit1_x,hit1_y,lay2_negx,lay2_negy);

      // center of positive curved track
      const float apos = 0.5*((hit0_x+hit1_x)+(hit0_y-hit1_y)*quad);
      const float bpos = 0.5*((hit0_y+hit1_y)-(hit0_x-hit1_x)*quad);
      
      // positive points of intersection with third layer
      float lay2_posx = 0.0f, lay2_posy = 0.0f;
      intersectThirdLayer(apos,bpos,hit1_x,hit1_y,lay2_posx,lay2_posy);

      const float lay2_z = 2*hit1_z-hit0_z; // for dz displacements -- straight line window
      const auto  lay2_etaBinM = getEtaPartition(getEta(lay2_r,lay2_z)-Config::dEtaSeedTrip);
      const auto  lay2_etaBinP = getEtaPartition(getEta(lay2_r,lay2_z)+Config::dEtaSeedTrip);
      const auto  lay2_phiBinM = getPhiPartition(getPhi(lay2_negx,lay2_negy));
      const auto  lay2_phiBinP = getPhiPartition(getPhi(lay2_posx,lay2_posy));
      
      std::vector<int> cand_hit2_indices = getCandHitIndices(lay2_etaBinM,lay2_etaBinP,lay2_phiBinM,lay2_phiBinP,segmentMap[2]);
      for (auto&& ihit2 : cand_hit2_indices){ // loop over candidate second layer hits
	const Hit & hit2 = evt_lay_hits[2][ihit2];

	const float lay1_predz = (hit0_z + hit2.z()) / 2.;
	// filter by residual of second layer hit
	if (std::abs(lay1_predz-hit1_z) > Config::seed_z1cut) continue;

	const float hit2_x = hit2.x();
	const float hit2_y = hit2.y();

	// now fit a circle, extract pT and d0 from center and radius
	const float mr = (hit1_y-hit0_y)/(hit1_x-hit0_x);
	const float mt = (hit2_y-hit1_y)/(hit2_x-hit1_x);
	const float a  = (mr*mt*(hit2_y-hit0_y) + mr*(hit1_x+hit2_x) - mt*(hit0_x+hit1_x))/(2.*(mr-mt));
	const float b  = -1.*(a-(hit0_x+hit1_x)/2.)/mr + (hit0_y+hit1_y)/2.;
	const float r  = getHypot(hit0_x-a,hit0_y-b);

	// filter by d0 cut 5mm, pT cut 0.5 GeV (radius of 0.5 GeV track)
	if ((r < Config::maxCurvR) || (std::abs(getHypot(a,b)-r) > Config::seed_d0cut)) continue; 
	
	// create a track object
	//int hitIndices[3] = {ihit0,ihit1,ihit2};
        // XXMT4K: Here will need layer indices, too now.
	// XX Track seed(dummystate,0.0f,seedID,Config::nlayers_per_seed,hitIndices); // argh! super not type-safe with dummystate

	// estimate and set the charge
	// XX seed.setCharge(calculateCharge(hit0,hit1,hit2));

	// save the track and track extra
	// XX evt_seed_tracks.push_back(seed);
	evt_seed_extras.emplace_back(seedID);

	// increment dummy counter for seedID
	seedID++; 
      } // end loop over third layer matches
    } // end loop over first layer matches
  } // end loop over second layer hits

  // now do conformal fit + KF fit
  fitSeeds(evt_lay_hits,evt_seed_tracks,ev);
}

void buildHitPairs(const std::vector<HitVec>& evt_lay_hits, const BinInfoLayerMap& segLayMap, PairIdxVec& hit_pairs){
  // use only average radius of inner radius for calculation
  // lay01angdiff is a parameter for phi search window derived numerically from Mathematica... see one of my old talks

  for (unsigned int ihit=0;ihit<evt_lay_hits[1].size();++ihit) { // 1 = second layer
    const float outerhitz = evt_lay_hits[1][ihit].z(); // remember, layer[0] is first layer! --> second layer = [1]
    const float outerphi  = evt_lay_hits[1][ihit].phi();

    const auto etaBinMinus = getEtaPartition(getEta(Config::fRadialSpacing,(outerhitz-Config::seed_z0cut)/2.));
    const auto etaBinPlus  = getEtaPartition(getEta(Config::fRadialSpacing,(outerhitz+Config::seed_z0cut)/2.));
    const auto phiBinMinus = getPhiPartition(normalizedPhi(outerphi - Config::lay01angdiff));
    const auto phiBinPlus  = getPhiPartition(normalizedPhi(outerphi + Config::lay01angdiff));

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

    const float quad = std::sqrt((4*maxCurvR2 - diffx2 -diffy2) / (diffx2 + diffy2));

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
    float posx2 = 0.0f, posy2 = 0;
    intersectThirdLayer(apos,bpos,x1,y1,posx2,posy2);
    const float posPhi = getPhi(posx2,posy2);

    // negative points of intersection with third layer
    float negx2 = 0.0f, negy2 = 0;
    intersectThirdLayer(aneg,bneg,x1,y1,negx2,negy2);
    const float negPhi = getPhi(negx2,negy2);

    const float thirdZline = 2*hit1.z()-hit0.z(); // for dz displacements -- straight line window
    const auto etaBinMinus = getEtaPartition(getEta(lay3rad,thirdZline)-Config::dEtaSeedTrip);
    const auto etaBinPlus  = getEtaPartition(getEta(lay3rad,thirdZline)+Config::dEtaSeedTrip);
    const auto phiBinMinus = getPhiPartition(negPhi);
    const auto phiBinPlus  = getPhiPartition(posPhi);

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

void intersectThirdLayer(const float a, const float b, const float hit1_x, const float hit1_y, float& lay2_x, float& lay2_y){
  const float a2 = a*a; const float b2 = b*b; const float a2b2 = a2+b2;
  const float lay2rad2  = (Config::fRadialSpacing*Config::fRadialSpacing)*9.0f; // average third radius squared
  const float maxCurvR2 = Config::maxCurvR * Config::maxCurvR;

  const float quad = std::sqrt( 2.0f*maxCurvR2*(a2b2+lay2rad2) - (a2b2-lay2rad2)*(a2b2-lay2rad2) - maxCurvR2*maxCurvR2 );
  const float pos[2] = { (a2*a + a*(b2+lay2rad2-maxCurvR2) - b*quad)/ a2b2 , (b2*b + b*(a2+lay2rad2-maxCurvR2) + a*quad)/ a2b2 };
  const float neg[2] = { (a2*a + a*(b2+lay2rad2-maxCurvR2) + b*quad)/ a2b2 , (b2*b + b*(a2+lay2rad2-maxCurvR2) - a*quad)/ a2b2 };

  // since we have two intersection points, arbitrate which one is closer to layer2 hit
  if (getHypot(pos[0]-hit1_x,pos[1]-hit1_y)<getHypot(neg[0]-hit1_x,neg[1]-hit1_y)) {
    lay2_x = pos[0];
    lay2_y = pos[1];
  }
  else {
    lay2_x = neg[0];
    lay2_y = neg[1];
  }
}

void buildHitTripletsApproxWindow(const std::vector<HitVec>& evt_lay_hits, const BinInfoLayerMap& segLayMap, 
				  const PairIdxVec& hit_pairs, TripletIdxVec& hit_triplets){
  const float thirdRad = Config::fRadialSpacing * 3.0; // average third radius

  for (auto&& hit_pair : hit_pairs){
    const Hit& hit0 = evt_lay_hits[0][hit_pair[0]];
    const Hit& hit1 = evt_lay_hits[1][hit_pair[1]];

    const float thirdZline = 2*hit1.z()-hit0.z(); // for dz displacements -- straight line window
    const auto etaBinMinus = getEtaPartition(getEta(thirdRad,thirdZline)-Config::dEtaSeedTrip);
    const auto etaBinPlus  = getEtaPartition(getEta(thirdRad,thirdZline)+Config::dEtaSeedTrip);

    const float linePhi = getPhi(hit1.position()[0] - hit0.position()[0], hit1.position()[1] - hit0.position()[1]);
    float thirdPhiMinus = 0.0f;
    float thirdPhiPlus  = 0.0f;  
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
    if ((r >= Config::maxCurvR) && (std::abs(getHypot(a,b)-r) <= Config::seed_d0cut)) {
      filtered_triplets.push_back(hit_triplet);
    } // d0 cut 5mm, pT cut 0.5 GeV (radius of 0.5 GeV track)
  }
}

void filterHitTripletsBySecondLayerZResidual(const std::vector<HitVec>& evt_lay_hits, const TripletIdxVec& hit_triplets, TripletIdxVec& filtered_triplets){
  for (auto&& hit_triplet : hit_triplets){
    //    const float z1 = predz(evt_lay_hits[0][hit_triplet[0]].z(),evt_lay_hits[0][hit_triplet[0]].r(),evt_lay_hits[2][hit_triplet[2]].z(),evt_lay_hits[2][hit_triplet[2]].r(),Config::fRadialSpacing*2.);
    const float z1 = (evt_lay_hits[0][hit_triplet[0]].z() + evt_lay_hits[2][hit_triplet[2]].z()) / 2.;
    if (std::abs(z1-evt_lay_hits[1][hit_triplet[1]].z()) < Config::seed_z1cut) { // three sigma inclusion
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
			    TrackVec& evt_seed_tracks, TrackExtraVec& evt_seed_extras, Event& ev)
{
  fprintf(stderr, "__FILE__::__LINE__ Needs fixing for B/E support, search for XXMT4K\n");
  exit(1);

  const PropagationFlags pflags(PF_none);

  // now perform kalman fit on seeds --> first need initial parameters --> get from Conformal fitter!
  unsigned int seedID = 0;
  for(auto&& hit_triplet : filtered_triplets){
    int charge = calculateCharge(evt_lay_hits[0][hit_triplet[0]],evt_lay_hits[1][hit_triplet[1]],evt_lay_hits[2][hit_triplet[2]]);

    TrackState updatedState;
    conformalFit(evt_lay_hits[0][hit_triplet[0]],evt_lay_hits[1][hit_triplet[1]],evt_lay_hits[2][hit_triplet[2]],updatedState,false); 
    updatedState.charge = charge;

    // CF is bad at getting a good pT estimate, phi and theta are fine
    // "best case" config found in other studies is to set TS by hand:
    // x,y,z to (0,0,0)
    // exx,eyy,ezz to 4*(Config::beamspotXY/Z)^2 (2 sigma)
    // px,py,pz from CF the same
    // epxpx,epypy,epzpz set to 0.25*(px/py/pz)^2 (half sigma)

    TSLayerPairVec updatedStates; // validation for position pulls
    float sumchi2 = 0;
    for (auto ilayer = 0; ilayer < Config::nlayers_per_seed; ++ilayer)
    {
      const Hit& seed_hit = evt_lay_hits[ilayer][hit_triplet[ilayer]];
      TrackState propState = propagateHelixToR(updatedState, seed_hit.r(), pflags);
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = seed_hit.measurementState();
      const float chi2 = computeChi2(propState,measState);
      sumchi2 += chi2;// --> could use this to make the chi2
      updatedState = updateParameters(propState, measState);
      updatedStates.push_back(std::make_pair(ilayer,updatedState)); // validation
    }

    ///int hitIndices[3] = {hit_triplet[0],hit_triplet[1],hit_triplet[2]};
    // XXMT4K: Here will need layer indices, too now.
    // XX Track seed(updatedState, sumchi2, seedID, Config::nlayers_per_seed, hitIndices);//fixme chi2
    // XX evt_seed_tracks.push_back(seed);
    evt_seed_extras.emplace_back(seedID);
    seedID++; // increment dummy counter for seedID
  }
}

void fitSeeds(const std::vector<HitVec>& evt_lay_hits, TrackVec& evt_seed_tracks, Event& ev){
  // copy+paste code to fit needs here...

  const PropagationFlags pflags(PF_none);

  for (auto&& seed : evt_seed_tracks)
  {
    // state to save to track
    TrackState updatedState;
    //const int seedID = seed.label(); // label == seedID!
    
    // first do CF Fit
    conformalFit(evt_lay_hits[0][seed.getHitIdx(0)],evt_lay_hits[1][seed.getHitIdx(1)],evt_lay_hits[2][seed.getHitIdx(2)],updatedState,false); 
    updatedState.charge = seed.charge();

    TSLayerPairVec updatedStates; // validation for position pulls
    float sumchi2 = 0.0f; // chi2 to save
    for (auto ilayer = 0; ilayer < Config::nlayers_per_seed; ++ilayer) {
      const Hit& seed_hit = evt_lay_hits[ilayer][seed.getHitIdx(ilayer)];
      TrackState propState = propagateHelixToR(updatedState, seed_hit.r(), pflags);
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = seed_hit.measurementState();
      const float chi2 = computeChi2(propState,measState);
      sumchi2 += chi2;// --> could use this to make the chi2
      updatedState = updateParameters(propState, measState);
      updatedStates.push_back(std::make_pair(ilayer,updatedState)); // validation
    }
    
    // make sure to set chi2 and track state of KF fit!
    seed.setChi2(sumchi2);
    seed.setState(updatedState);
  }
}
