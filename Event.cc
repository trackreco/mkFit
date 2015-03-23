#include "Event.h"
#include "Simulation.h"
#include "KalmanUtils.h"
//#include "seedtest.h"
#include "buildtest.h"
#include "fittest.h"
#include "ConformalUtils.h"
#include "Debug.h"

#ifdef TBB
#include "tbb/tbb.h"
#endif

static bool sortByPhi(const Hit& hit1, const Hit& hit2)
{
  return hit1.phi()<hit2.phi();
}

static bool tracksByPhi(const Track& t1, const Track& t2)
{
  return t1.posPhi()<t2.posPhi();
}

#ifdef ETASEG
const float etaDet = 2.0;
static bool sortByEta(const Hit& hit1, const Hit& hit2){
  return hit1.eta()<hit2.eta();
}

// within a layer with a "reasonable" geometry, ordering by Z is the same as eta
static bool sortByZ(const Hit& hit1, const Hit& hit2){
  return hit1.z()<hit2.z();
}
#endif

Event::Event(const Geometry& g, Validation& v, int threads) : geom_(g), validation_(v), threads_(threads)
{
  layerHits_.resize(geom_.CountLayers());
  segmentMap_.resize(geom_.CountLayers());
}

void Event::Simulate(unsigned int nTracks)
{
  simTracks_.resize(nTracks);
  for (auto&& l : layerHits_) {
    l.reserve(nTracks);
  }

#ifdef TBB
  parallel_for( tbb::blocked_range<size_t>(0, nTracks, 100), 
      [&](const tbb::blocked_range<size_t>& itracks) {

    const Geometry tmpgeom(geom_.clone()); // USolids isn't thread safe??
    for (auto itrack = itracks.begin(); itrack != itracks.end(); ++itrack) {
#else
    const Geometry& tmpgeom(geom_);
    for (unsigned int itrack=0; itrack<nTracks; ++itrack) {
#endif
      //create the simulated track
      SVector3 pos;
      SVector3 mom;
      SMatrixSym66 covtrk;
      HitVec hits, initialhits;
      // unsigned int starting_layer  = 0; --> for displaced tracks, may want to consider running a separate Simulate() block with extra parameters

      int q=0;//set it in setup function
      float pt = 0.5+g_unif(g_gen)*9.5;//this input, 0.5<pt<10 GeV (below ~0.5 GeV does not make 10 layers)
      setupTrackByToyMC(pos,mom,covtrk,hits,itrack,q,pt,tmpgeom,initialhits);
      Track sim_track(q,pos,mom,covtrk,hits,0,initialhits);
      simTracks_[itrack] = sim_track;
    }
#ifdef TBB
  });
#endif

  // fill vector of hits in each layer
  for (const auto& track : simTracks_) {
    for (const auto& hit : track.hitsVector()) {
      layerHits_[hit.layer()].push_back(hit);
    }
  }
  validation_.fillSimHists(simTracks_);
}

void Event::Segment()
{
#ifdef DEBUG
  bool debug=true;
#endif
  //sort in phi and dump hits per layer, fill phi partitioning
  for (unsigned int ilayer=0; ilayer<layerHits_.size(); ++ilayer) {
    dprint("Hits in layer=" << ilayer);
    
#ifdef ETASEG
    segmentMap_[ilayer].resize(Config::nEtaPart);    
    // eta first then phi
    std::sort(layerHits_[ilayer].begin(), layerHits_[ilayer].end(), sortByZ);
    std::vector<unsigned int> lay_eta_bin_count(Config::nEtaPart);
    for (unsigned int ihit=0;ihit<layerHits_[ilayer].size();++ihit) {
      unsigned int etabin = getEtaPartition(layerHits_[ilayer][ihit].eta(),etaDet);
      dprint("ihit: " << ihit << " eta: " << layerHits_[ilayer][ihit].eta() << " etabin: " << etabin);
      lay_eta_bin_count[etabin]++;
    }
    //now set index and size in partitioning map and then sort the bin by phi
    
    int lastEtaIdxFound = -1;
    int lastPhiIdxFound = -1;

    for (unsigned int etabin=0; etabin<Config::nEtaPart; ++etabin) {
      unsigned int firstEtaBinIdx = lastEtaIdxFound+1;
      unsigned int etaBinSize = lay_eta_bin_count[etabin];
      if (etaBinSize>0){
        lastEtaIdxFound+=etaBinSize;
      }

      //sort by phi in each "eta bin"
      std::sort(layerHits_[ilayer].begin() + firstEtaBinIdx,layerHits_[ilayer].begin() + (etaBinSize+firstEtaBinIdx), sortByPhi); // sort from first to last in eta
      std::vector<unsigned int> lay_eta_phi_bin_count(Config::nPhiPart);

      for(unsigned int ihit = firstEtaBinIdx; ihit < etaBinSize+firstEtaBinIdx; ++ihit){
        dprint("ihit: " << ihit << " r(layer): " << layerHits_[ilayer][ihit].r() << "(" << ilayer << ") phi: " 
                        << layerHits_[ilayer][ihit].phi() << " phipart: " << getPhiPartition(layerHits_[ilayer][ihit].phi()) << " eta: "
                        << layerHits_[ilayer][ihit].eta() << " etapart: " << getEtaPartition(layerHits_[ilayer][ihit].eta(),etaDet));
        unsigned int phibin = getPhiPartition(layerHits_[ilayer][ihit].phi());
        lay_eta_phi_bin_count[phibin]++;
      }

      for (unsigned int phibin=0; phibin<Config::nPhiPart; ++phibin) {
        unsigned int firstPhiBinIdx = lastPhiIdxFound+1;
        unsigned int phiBinSize = lay_eta_phi_bin_count[phibin];
        BinInfo phiBinInfo(firstPhiBinIdx,phiBinSize);
        segmentMap_[ilayer][etabin].push_back(phiBinInfo);
        if (phiBinSize>0){
          lastPhiIdxFound+=phiBinSize;
        }
#ifdef DEBUG
        if ((debug) && (phiBinSize !=0)) printf("ilayer: %1u etabin: %1u phibin: %2u first: %2u last: %2u \n", 
                                                ilayer, etabin, phibin, 
                                                segmentMap_[ilayer][etabin][phibin].first, 
                                                segmentMap_[ilayer][etabin][phibin].second+segmentMap_[ilayer][etabin][phibin].first
                                                );
#endif
      } // end loop over storing phi index
    } // end loop over storing eta index
#else
    segmentMap_[ilayer].resize(1);    // only one eta bin for special case, avoid ifdefs
    std::sort(layerHits_[ilayer].begin(), layerHits_[ilayer].end(), sortByPhi);
    std::vector<unsigned int> lay_phi_bin_count(Config::nPhiPart);//should it be 63? - yes!
    for (unsigned int ihit=0;ihit<layerHits_[ilayer].size();++ihit) {
      dprint("hit r/phi/eta : " << layerHits_[ilayer][ihit].r() << " "
                                << layerHits_[ilayer][ihit].phi() << " " << layerHits_[ilayer][ihit].eta());

      unsigned int phibin = getPhiPartition(layerHits_[ilayer][ihit].phi());
      lay_phi_bin_count[phibin]++;
    }

    //now set index and size in partitioning map
    int lastIdxFound = -1;
    for (unsigned int bin=0; bin<Config::nPhiPart; ++bin) {
      unsigned int binSize = lay_phi_bin_count[bin];
      unsigned int firstBinIdx = lastIdxFound+1;
      BinInfo binInfo(firstBinIdx, binSize);
      segmentMap_[ilayer][0].push_back(binInfo); // [0] bin is just the only eta bin ... reduce ifdefs
      if (binSize>0){
        lastIdxFound+=binSize;
      }
    }
#endif
  } // end loop over layers
}

void Event::Seed()
{
#define SIMSEEDS
#ifdef SIMSEEDS
  //create seeds (from sim tracks for now)
  for (unsigned int itrack=0;itrack<simTracks_.size();++itrack) {
    const Track& trk = simTracks_[itrack];
    const HitVec& hits = trk.hitsVector();
    HitVec seedhits;
    float chi2 = 0;
    
    //#define CONF_SEED
#ifdef CONF_SEED
    int charge = trk.charge();

    //    std::cout << "charge: " << charge << std::endl;
    //    if (hits[1].phi() > hits[2].phi()){charge = 1;}
    //    else {charge = -1;}

    std::cout << "Conformal track: " << itrack << std::endl;

    TrackState cfitStateHit0;
    conformalFit(hits[0],hits[1],hits[2],charge,cfitStateHit0);

    TrackState updatedState = cfitStateHit0;
    updatedState.errors*=10.;
#else
    TrackState updatedState = trk.state();
#endif // CONFSEED
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
      seedhits.push_back(seed_hit);//fixme chi2
      chi2 += computeChi2(propState,measState);
    }
    Track seed(updatedState,seedhits,chi2);//fixme chi2
    seedTracks_.push_back(seed);
  }
#else
  //  runSeedTest(*this);
  bool debug = false;

  // follow CMSSW -> start with hits in 2nd layer, build seed by then projecting back to beamspot.
  // build seed pairs... then project to third layer.

  // p=qBR => R is rad of curvature (in meters), q = 0.2997 in natural units, B = 3.8T, p = min pt = 0.5 GeV. 
  // so R = (0.5/0.2997...*3.8) * 100 -> max rad. of curv. in cm 

  std::vector<HitVec> seed_pairs;

  //  const float curvRad = 43.8900125; // 0.5 GeV track
  //  const float d0 = 0.1; // 1 mm x,y beamspot from simulation
  const float dZ = 1.0; // 3.25 - 3.3
  const float innerrad = 4.0; // average radius of inner radius

  const float alphaBeta = 0.0520195; // 0.0458378 --> for d0 = .0025 cm
  // min alphaBeta for 100% eff Pi/32 - 0.02  to Pi/2 - 0.025 

  // 0.5 GeV track d0 = 0   has eff = 48410 /50000, fakes = 203.2k, aB = 0.0456793
  // 0.5 GeV track d0 = 0.1 has eff = 48422 /50000, fakes = 217.1k, aB = 0.0520195

  // 10  GeV track d0 = 0   has eff = 43856 /50000, fakes = 110.8k, aB = 0.00227844
  // 10  GeV Track d0 = 0.1 has eff = 45915 /50000, fakes = 124.1k, aB = 0.00852886

  // Full acceptance for efficiency needs alphaBeta =Pi/32 - 0.2 &&  dZ = 3.3.  Fakes at 472k
  // At nominal parameters alphaBeta = 0.0520195 && dZ = 1.0, eff = 48422 / 50000, fakes at 217k

  for (unsigned int ihit=0;ihit<layerHits_[1].size();++ihit) { // 1 = second layer
    float outerrad  = layerHits_[1][ihit].r();
    float outerphi  = layerHits_[1][ihit].phi();
    float outerhitz = layerHits_[1][ihit].position()[2];
    
    // for d0 displacements -- helix phi window

    //    float alphaPlus  = acos(-((d0*d0)+(outerrad*outerrad)-2.*(d0*curvRad))/((2.*outerrad)*(d0-curvRad)));
    //   float betaPlus   = acos(-((d0*d0)+(innerrad*innerrad)-2.*(d0*curvRad))/((2.*innerrad)*(d0-curvRad)));

    //    float alphaMinus = acos(-((d0*d0)+(outerrad*outerrad)-2.*(d0*curvRad))/((2.*outerrad)*(-d0+curvRad)));
    //    float betaMinus  = acos(-((d0*d0)+(innerrad*innerrad)-2.*(d0*curvRad))/((2.*innerrad)*(-d0+curvRad)));

    //    float innerPhiMinus = normalizedPhi(outerphi - alphaMinus + betaMinus);
    //   float innerPhiPlus  = normalizedPhi(outerphi - alphaPlus  + betaPlus);

    float innerPhiMinus = normalizedPhi(outerphi - alphaBeta);
    float innerPhiPlus  = normalizedPhi(outerphi + alphaBeta);

    unsigned int phiBinMinus = getPhiPartition(innerPhiMinus);
    unsigned int phiBinPlus  = getPhiPartition(innerPhiPlus);

    // for dz displacements -- straight line window
    
    //    float innerZPlus  = (innerrad/outerrad)*(outerhitz-dZ)+dZ;
    //    float innerZMinus = (innerrad/outerrad)*(outerhitz+dZ)-dZ;
    float innerZPlus  = (0.5)*(outerhitz-dZ)+dZ;
    float innerZMinus = (0.5)*(outerhitz+dZ)-dZ;
    
#ifdef ETASEG
    unsigned int etaBinMinus = getEtaPartition(normalizedEta(getEta(innerrad,innerZMinus)),etaDet);
    unsigned int etaBinPlus  = getEtaPartition(normalizedEta(getEta(innerrad,innerZPlus)),etaDet);
#else
    unsigned int etaBinMinus = 0;
    unsigned int etaBinPlus  = 0;
#endif    

    unsigned int firstPhiPart = getPhiPartition(simTracks_[layerHits_[1][ihit].mcTrackID()].hitsVector()[0].phi());
    unsigned int firstEtaPart = getEtaPartition(simTracks_[layerHits_[1][ihit].mcTrackID()].hitsVector()[0].eta(),etaDet);
    /*    if (phiBinMinus<=phiBinPlus){
      if ( (firstPhiPart < phiBinMinus) || (firstPhiPart > phiBinPlus) ){
	missed_phis.push_back(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2]);
      }
    }
    else{
      if ( (firstPhiPart > phiBinMinus) || (firstPhiPart < phiBinPlus) ){
	missed_phis.push_back(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2]);
      }
    }
    if ( (firstEtaPart < etaBinMinus) || (firstEtaPart > etaBinPlus) ){
      missed_etas.push_back(simTracks_[seed_pair[0].mcTrackID()].hitsVector()[2]);
    }
    
    */

    if (debug) {
      unsigned int mcID = layerHits_[1][ihit].mcTrackID();
      Hit innerHit = simTracks_[mcID].hitsVector()[0];      
      float innerhitz = innerHit.position()[2];
      float innerphi  = innerHit.phi();
      printf("ihit: %1u iphi: %08.5f \n   phiM: %08.5f phiC: %08.5f phiP: %08.5f \n   binM: %2u binC: %2u binP: %2u \n",
	     ihit,innerphi, 
	     innerPhiMinus,outerphi,innerPhiPlus,
	     phiBinMinus,getPhiPartition(outerphi),phiBinPlus
	     );
  
      printf("   innerZ: %010.5f \n   ZM: %010.5f ZC: %010.5f ZP: %010.5f \n",
	     innerhitz,
	     innerZMinus,(innerrad/outerrad)*outerhitz,innerZPlus
	     );
    }

    for (unsigned int ieta = etaBinMinus; ieta <= etaBinPlus; ++ieta){
      
      BinInfo binInfoMinus = lay_eta_phi_hit_idx_[0][ieta][int(phiBinMinus)];
      BinInfo binInfoPlus  = lay_eta_phi_hit_idx_[0][ieta][int(phiBinPlus)];
      
      std::vector<unsigned int> cand_hit_idx;
      std::vector<unsigned int>::iterator index_iter; // iterator for vector
      
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
	unsigned int etaBinSize = lay_eta_phi_hit_idx_[0][ieta][62].first+lay_eta_phi_hit_idx_[0][ieta][62].second;
	
	for (unsigned int ihit  = firstIndex; ihit < etaBinSize; ++ihit){
	  cand_hit_idx.push_back(ihit);
	}
	
	unsigned int etaBinStart= lay_eta_phi_hit_idx_[0][ieta][0].first;
	unsigned int maxIndex   = binInfoPlus.first+binInfoPlus.second;
	
	for (unsigned int ihit  = etaBinStart; ihit < maxIndex; ++ihit){
	  cand_hit_idx.push_back(ihit);
	}
      }

      // build a pair for each 2nd layer hit
      
      for (index_iter = cand_hit_idx.begin(); index_iter != cand_hit_idx.end(); ++index_iter){
	HitVec seed_pair;
	seed_pair.push_back(layerHits_[0][*index_iter]);
	seed_pair.push_back(layerHits_[1][ihit]);
	seed_pairs.push_back(seed_pair);
      }
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
      updatedState = updateParameters(propState, measState,projMatrix36_,projMatrix36T_);
    }
    Track seed(updatedState,seed_triplet,0.);//fixme chi2
    seedTracks_.push_back(seed);
  }

  TrackVec seedTracksMC;

  for (unsigned int itrack=0;itrack<simTracks_.size();++itrack) {
    HitVec seedhits;
    Track& trk = simTracks_[itrack];
    HitVec& hits = trk.hitsVector();
    float chi2   = 0;
    TrackState updatedState = trk.state();

    for (auto ilayer=0U;ilayer<nlayers_per_seed;++ilayer) {//seeds have first three layers as seeds
      Hit seed_hit = hits[ilayer]; // do this for now to make initHits element number line up with HitId number
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
#ifdef CHECKSTATEVALID
      if (!propState.valid) {
        break;
      }
#endif
      MeasurementState measState = seed_hit.measurementState();
      updatedState = updateParameters(propState,measState);
      seedhits.push_back(seed_hit);//fixme chi2
      chi2 += computeChi2(propState,measState);
    }
    Track seed(updatedState,seedhits,chi2);//fixme chi2
    seedTracksMC.push_back(seed);
  }

  validation_.fillSeedHists(seed_pairs,seed_triplets,simTracks_,missed_phis,missed_etas,deta,chi2triplet,conformalTracks,seedTracks_,seedTracksMC);
#endif
  std::sort(seedTracks_.begin(), seedTracks_.end(), tracksByPhi);
}

void Event::Find()
{
  //buildTracksBySeeds(*this);
  buildTracksByLayers(*this);
  validation_.fillAssociationHists(candidateTracks_,simTracks_);
  validation_.fillCandidateHists(candidateTracks_);
}

void Event::Fit()
{
#ifdef ENDTOEND
  runFittingTest(*this, candidateTracks_);
#else
  runFittingTest(*this, simTracks_);
#endif
}

