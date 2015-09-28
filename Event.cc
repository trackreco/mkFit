#include "Event.h"
#include "Simulation.h"
#include "KalmanUtils.h"
#include "seedtest.h"
#include "buildtest.h"
#include "fittest.h"
#include "BinInfoUtils.h"
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
static bool sortByEta(const Hit& hit1, const Hit& hit2){
  return hit1.eta()<hit2.eta();
}

// within a layer with a "reasonable" geometry, ordering by Z is the same as eta
static bool sortByZ(const Hit& hit1, const Hit& hit2){
  return hit1.z()<hit2.z();
}
#endif

Event::Event(const Geometry& g, Validation& v, unsigned int evtID, int threads) : geom_(g), validation_(v), evtID_(evtID), threads_(threads)
{
  layerHits_.resize(Config::nLayers);
  segmentMap_.resize(Config::nLayers);

  validation_.resetValidationMaps(); // need to reset maps for every event.
}

void Event::Simulate()
{
  simTracks_.resize(Config::nTracks);
  for (auto&& l : layerHits_) {
    l.reserve(Config::nTracks);
  }

#ifdef TBB
  parallel_for( tbb::blocked_range<size_t>(0, Config::nTracks, 100), 
      [&](const tbb::blocked_range<size_t>& itracks) {

    const Geometry tmpgeom(geom_.clone()); // USolids isn't thread safe??
    for (auto itrack = itracks.begin(); itrack != itracks.end(); ++itrack) {
#else
    const Geometry& tmpgeom(geom_);
    for (unsigned int itrack=0; itrack<Config::nTracks; ++itrack) {
#endif
      //create the simulated track
      SVector3 pos;
      SVector3 mom;
      SMatrixSym66 covtrk;
      HitVec hits;
      TSVec  initialTSs;
      // unsigned int starting_layer  = 0; --> for displaced tracks, may want to consider running a separate Simulate() block with extra parameters

      int q=0;//set it in setup function
      setupTrackByToyMC(pos,mom,covtrk,hits,itrack,q,tmpgeom,initialTSs); // do the simulation
      validation_.collectSimTkTSVecMapInfo(itrack,initialTSs); // save initial TS parameters
      Track sim_track(q,pos,mom,covtrk,hits,0,itrack,itrack);
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
      unsigned int etabin = getEtaPartition(layerHits_[ilayer][ihit].eta());
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
	                << layerHits_[ilayer][ihit].eta() << " etapart: " << getEtaPartition(layerHits_[ilayer][ihit].eta()));
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

#ifdef DEBUG
  for (unsigned int ilayer = 0; ilayer < Config::nLayers; ilayer++) {
    unsigned int etahitstotal = 0;
    for (unsigned int etabin = 0; etabin < Config::nEtaPart; etabin++){
      unsigned int etahits = segmentMap_[ilayer][etabin][Config::nPhiPart-1].first + segmentMap_[ilayer][etabin][Config::nPhiPart-1].second - segmentMap_[ilayer][etabin][0].first;
      std::cout << "etabin: " << etabin << " hits in bin: " << etahits << std::endl;
      etahitstotal += etahits;

      for (unsigned int phibin = 0; phibin < Config::nPhiPart; phibin++){
	//	if (segmentMap_[ilayer][etabin][phibin].second > 3) {std::cout << "   phibin: " << phibin << " hits: " << segmentMap_[ilayer][etabin][phibin].second << std::endl;}
      }
    }
    std::cout << "layer: " << ilayer << " totalhits: " << etahitstotal << std::endl;
  }
#endif
}

void Event::Seed()
{
#ifdef ENDTOEND
  //  buildSeedsByRoadTriplets(layerHits_,segmentMap_,seedTracks_,*this);
  buildSeedsByMC(simTracks_,seedTracks_,*this);
#else
  buildSeedsByMC(simTracks_,seedTracks_,*this);
#endif
  std::sort(seedTracks_.begin(), seedTracks_.end(), tracksByPhi);
}

void Event::Find()
{
  buildTracksBySeeds(*this);
  //  buildTracksByLayers(*this);

  // From CHEP-2015
  // buildTestSerial(*this, Config::nlayers_per_seed, Config::maxCand, Config::chi2Cut, Config::nSigma, Config::minDPhi);
}

void Event::Fit()
{
#ifdef ENDTOEND
  runFittingTest(*this, candidateTracks_);
#else
  runFittingTest(*this, simTracks_);
#endif
}

void Event::Validate(const unsigned int ievt){
  validation_.fillSegmentTree(segmentMap_,ievt);
  validation_.fillBranchTree(ievt);
  validation_.makeSimTkToRecoTksMaps(seedTracks_,candidateTracks_,fitTracks_);
  validation_.fillEffTree(simTracks_,ievt);
  validation_.makeSeedTkToRecoTkMaps(candidateTracks_,fitTracks_);
  validation_.fillFakeRateTree(seedTracks_,ievt);
}
