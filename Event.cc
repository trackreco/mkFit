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

inline bool sortByPhi(const Hit& hit1, const Hit& hit2)
{
  return hit1.phi()<hit2.phi();
}

static bool tracksByPhi(const Track& t1, const Track& t2)
{
  return t1.posPhi()<t2.posPhi();
}

#ifdef ETASEG
inline bool sortByEta(const Hit& hit1, const Hit& hit2){
  return hit1.eta()<hit2.eta();
}

// within a layer with a "reasonable" geometry, ordering by Z is the same as eta
inline bool sortByZ(const Hit& hit1, const Hit& hit2){
  return hit1.z()<hit2.z();
}
#endif

void Event::resetLayerHitMap(bool resetSimHits) {
  layerHitMap_.clear();
  layerHitMap_.resize(simHitsInfo_.size());
  for (int ilayer = 0; ilayer < layerHits_.size(); ++ilayer) {
    for (int index = 0; index < layerHits_[ilayer].size(); ++index) {
      auto& hit = layerHits_[ilayer][index];
      layerHitMap_[hit.mcHitID()] = HitID(ilayer, index);
    }
  }
  if (resetSimHits) {
    for (auto&& track : simTracks_) {
      for (int il = 0; il < track.nTotalHits(); ++il) {
        track.setHitIdx(il, layerHitMap_[track.getHitIdx(il)].index);
      }
    }
  }
}

Event::Event(const Geometry& g, Validation& v, unsigned int evtID, int threads) : geom_(g), validation_(v), evtID_(evtID), threads_(threads)
{
  layerHits_.resize(Config::nLayers);
  segmentMap_.resize(Config::nLayers);

  validation_.resetValidationMaps(); // need to reset maps for every event.
}

void Event::Simulate()
{
  MCHitInfo::mcHitIDCounter_ = 0;

  simTracks_.resize(Config::nTracks);
  simHitsInfo_.resize(Config::nTotHit * Config::nTracks);
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
      setupTrackByToyMC(pos,mom,covtrk,hits,simHitsInfo_,itrack,q,tmpgeom,initialTSs); // do the simulation
      //setupTrackFromTextFile(pos,mom,covtrk,hits,simHitsInfo_,itrack,q,tmpgeom,initialTSs);
      validation_.collectSimTkTSVecMapInfo(itrack,initialTSs); // save initial TS parameters

      simTracks_[itrack] = Track(q,pos,mom,covtrk,0.0f);
      auto& sim_track = simTracks_[itrack];
      sim_track.setLabel(itrack);
      for (int ilay = 0; ilay < hits.size(); ++ilay) {
        sim_track.addHitIdx(hits[ilay].mcHitID(),0.0f); // tmp because of sorting...
        layerHits_[ilay].push_back(hits[ilay]);
      }
    }
#ifdef TBB
  });
#endif
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
  resetLayerHitMap(true);
}

void Event::Seed()
{
#ifdef ENDTOEND
  //  buildSeedsByRoadTriplets(layerHits_,segmentMap_,seedTracks_,*this);
  buildSeedsByMC(simTracks_,seedTracks_,seedTracksExtra_,*this);
#else
  buildSeedsByMC(simTracks_,seedTracks_,seedTracksExtra_,*this);
#endif
  // if we sort here, also have to sort seedTracksExtra and redo labels.
  //std::sort(seedTracks_.begin(), seedTracks_.end(), tracksByPhi);
}

void Event::Find()
{
  buildTracksBySeeds(*this);
  //  buildTracksByLayers(*this);

  // From CHEP-2015
  // buildTestSerial(*this, Config::nlayers_per_seed, Config::maxCandsPerSeed, Config::chi2Cut, Config::nSigma, Config::minDPhi);
}

void Event::Fit()
{
  fitTracks_.resize(candidateTracks_.size());
  fitTracksExtra_.resize(candidateTracks_.size());
#ifdef ENDTOEND
  runFittingTest(*this, candidateTracks_, candidateTracksExtra_);
#else
  runFittingTest(*this, simTracks_, simTracksExtra_);
#endif
}

void Event::Validate(const unsigned int ievt){
  validation_.fillSegmentTree(segmentMap_,ievt);
  validation_.fillBranchTree(ievt);
  validation_.makeSimTkToRecoTksMaps(*this);
  validation_.fillEffTree(simTracks_,ievt);
  validation_.makeSeedTkToRecoTkMaps(candidateTracks_,fitTracks_);
  validation_.fillFakeRateTree(seedTracks_,ievt);
}

void Event::PrintStats(const TrackVec& trks, TrackExtraVec& trkextras)
{
  int miss(0), found(0), fp_10(0), fp_20(0), hit8(0), h8_10(0), h8_20(0);

  for (auto&& trk : trks) {
    auto&& extra = trkextras[trk.label()];
    extra.setMCTrackIDInfo(trk, layerHits_, simHitsInfo_);
    if (extra.isMissed()) {
      ++miss;
    } else {
      auto&& mctrk = simTracks_[extra.mcTrackID()];
      auto pr = trk.pT()/mctrk.pT();
      found++;
      bool h8 = trk.nFoundHits() >= 8;
      bool pt10 = pr > 0.9 && pr < 1.1;
      bool pt20 = pr > 0.8 && pr < 1.2;
      fp_10 += pt10;
      fp_20 += pt20;
      hit8 += h8;
      h8_10 += h8 && pt10;
      h8_20 += h8 && pt20;
    }
  }
  std::cout << "found tracks=" << found   << "  in pT 10%=" << fp_10    << "  in pT 20%=" << fp_20    << "     no_mc_assoc="<< miss <<std::endl
            << "  nH >= 8   =" << hit8    << "  in pT 10%=" << h8_10    << "  in pT 20%=" << h8_20    << std::endl;
}
