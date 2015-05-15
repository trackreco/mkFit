#include "Event.h"
#include "Simulation.h"
#include "KalmanUtils.h"
#include "buildtest.h"
#include "fittest.h"
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
    TrackState updatedState = trk.state();
    HitVec seedhits;

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
    }
    Track seed(updatedState,seedhits,0.);//fixme chi2
    seedTracks_.push_back(seed);
  }
#else

  // follow CMSSW -> start with hits in 2nd layer, build seed by then projecting back to beamspot.
  // build seed pairs... then project to third layer.

  // p=qBR => R is rad of curvature (in meters), q = 0.2997 in natural units, B = 3.8T, p = min pt = 0.5 GeV. 
  // so R = (0.5/0.2997...*3.8) * 100 -> max rad. of curv. in cm 

  const float curvRad = 4.38900125;
  const float d0 = 0.1; // 1 mm x,y beamspot from simulation
  const float dZ = 2.0;

  for (unsigned int ihit=0;ihit<layerHits_[1].size();++ihit) { // 1 = second layer
    float outerrad  = layerHits_[1][ihit].r();
    float ouerphi   = layerHits_[1][ihit].phi();

    unsigned int mcID = layerHits_[1][ihit].mcTrackID();
    Hit innerHit = simTracks_[mcID].hitsVector()[0];
    float innerrad  = innerHit.r();
    innerrad = 4.0;
  
    std::cout << "Diff: " << innerrad - 4.0 <<std::endl;

    float innerPhiPlusCentral  = outerphi-acos(outerrad/(2.*curvRad))+acos(innerrad/(2.*curvRad));
    float innerPhiMinusCentral = outerphi-acos(outerrad/(-2.*curvRad))+acos(-innerrad/(2.*curvRad));

    // for d0 displacements

    float alphaPlus = acos(-((d0*d0)+(outerrad*outerrad)-2.*(d0*curvRad))/((2.*outerrad)*(d0-curvRad)));
    float betaPlus  = acos(-((d0*d0)+(innerrad*innerrad)-2.*(d0*curvRad))/((2.*innerrad)*(d0-curvRad)));

    float alphaMinus = acos(-((d0*d0)+(outerrad*outerrad)+2.*(d0*curvRad))/((2.*outerrad)*(d0+curvRad)));
    float betaMinus  = acos(-((d0*d0)+(innerrad*innerrad)+2.*(d0*curvRad))/((2.*innerrad)*(d0+curvRad)));

    float innerPhiPlus = outerphi-alphaPlus+betaPlus;
    float innerPhiMinus = outerphi-alphaMinus+betaMinus;

    float innerZPlus  = (innerrad/outerrad)*(outerhitz-dZ)+dZ;
    float innerZMinus = (innerrad/outerrad)*(outerhitz+dZ)-dZ;
    float centralZ    = (innerrad/outerrad)*outerhitz;

    printf("ihit: %1u \n   iphi: %5f ophi: %5f \n   iphiP: %5f iphiM: %5f \n   iphiPC: %5f iphiMC: %5f \n",
           ihit,
           innerHit.phi(), outerphi,
           innerPhiPlus,innerPhiMinus,
           innerPhiPlusCentral,innerPhiMinusCentral
           );
    printf("   innerZ: %5f iZM: %5f iZP: %5f cZ: %5f \n",
           innerhitz,innerZMinus,innerZPlus,centralZ
           );

    unsigned int phibinPlus  = getPhiPartition(innerPhiPlus);
    unsigned int phibinMinus = getPhiPartition(innerPhiMinus);

  }
#endif
  std::sort(seedTracks_.begin(), seedTracks_.end(), tracksByPhi);
}

void Event::Find()
{
  //buildTracksBySeeds(*this);
  buildTracksByLayers(*this);

  // From CHEP-2015
  // buildTestSerial(*this, Config::nlayers_per_seed, Config::maxCand, Config::chi2Cut, Config::nSigma, Config::minDPhi);

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

