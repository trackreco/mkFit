#include "Event.h"
#include "Simulation.h"
#include "KalmanUtils.h"
#include "seedtest.h"
#include "buildtest.h"
#include "fittest.h"
#include "BinInfoUtils.h"
#include "ConformalUtils.h"

//#define DEBUG
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
  //gc: not sure what is being done here
  layerHitMap_.clear();
  layerHitMap_.resize(simHitsInfo_.size());
  for (int ilayer = 0; ilayer < layerHits_.size(); ++ilayer) {
    for (int index = 0; index < layerHits_[ilayer].size(); ++index) {
      auto& hit = layerHits_[ilayer][index];
      assert(hit.mcHitID() >= 0); // tmp debug
      assert(hit.mcHitID() < layerHitMap_.size());
      layerHitMap_[hit.mcHitID()] = HitID(ilayer, index);
    }
  }
  if (resetSimHits) {
    for (auto&& track : simTracks_) {
      for (int il = 0; il < track.nTotalHits(); ++il) {
        assert(layerHitMap_[track.getHitIdx(il)].index >= 0); // tmp debug
        track.setHitIdx(il, layerHitMap_[track.getHitIdx(il)].index);
      }
    }
  }
}

Event::Event(const Geometry& g, Validation& v, int evtID, int threads) : geom_(g), validation_(v), evtID_(evtID), threads_(threads)
{
  layerHits_.resize(Config::nLayers);
  segmentMap_.resize(Config::nLayers);

  validation_.resetValidationMaps(); // need to reset maps for every event.
  if (Config::super_debug) {
    validation_.resetDebugVectors(); // need to reset vectors for every event.
  }
}

void Event::Simulate()
{
  MCHitInfo::mcHitIDCounter_ = 0;

  simTracks_.resize(Config::nTracks);
  simHitsInfo_.resize(Config::nTotHit * Config::nTracks);
  for (auto&& l : layerHits_) {
    l.resize(Config::nTracks);  // thread safety
  }
  simTrackStates_.resize(Config::nTracks);

#ifdef TBB
  parallel_for( tbb::blocked_range<size_t>(0, Config::nTracks, 100), 
      [&](const tbb::blocked_range<size_t>& itracks) {

    const Geometry tmpgeom(geom_.clone()); // USolids isn't thread safe??
    for (auto itrack = itracks.begin(); itrack != itracks.end(); ++itrack) {
#else
    const Geometry& tmpgeom(geom_);
    for (int itrack=0; itrack<Config::nTracks; ++itrack) {
#endif
      //create the simulated track
      SVector3 pos;
      SVector3 mom;
      SMatrixSym66 covtrk;
      HitVec hits;
      TSVec  initialTSs;
      // int starting_layer  = 0; --> for displaced tracks, may want to consider running a separate Simulate() block with extra parameters

      int q=0;//set it in setup function
      // do the simulation
      if (Config::useCMSGeom) setupTrackFromTextFile(pos,mom,covtrk,hits,simHitsInfo_,itrack,q,tmpgeom,initialTSs);
      else if (Config::endcapTest) setupTrackByToyMCEndcap(pos,mom,covtrk,hits,simHitsInfo_,itrack,q,tmpgeom,initialTSs);
      else setupTrackByToyMC(pos,mom,covtrk,hits,simHitsInfo_,itrack,q,tmpgeom,initialTSs); 

#ifdef CCSCOORD
      // setupTrackByToyMCEndcap is already in CCS coord, no need to convert
      if (Config::endcapTest==false) {
	float pt = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
	mom=SVector3(1./pt,atan2(mom[1],mom[0]),atan2(pt,mom[2]));
	for (int its = 0; its < initialTSs.size(); its++){
	  initialTSs[its].convertFromCartesianToCCS();
	}
      }
#endif
      // uber ugly way of getting around read-in / write-out of objects needed for validation
      if (Config::normal_val || Config::fit_val) {simTrackStates_[itrack] = initialTSs;}
      validation_.collectSimTkTSVecMapInfo(itrack,initialTSs); // save initial TS parameters in validation object ... just a copy of the above line

      simTracks_[itrack] = Track(q,pos,mom,covtrk,0.0f);
      auto& sim_track = simTracks_[itrack];
      sim_track.setLabel(itrack);
      for (int ilay = 0; ilay < hits.size(); ++ilay) {
        sim_track.addHitIdx(hits[ilay].mcHitID(),0.0f); // set to the correct hit index after sorting
        layerHits_[ilay][itrack] = hits[ilay];  // thread safety
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
  for (int ilayer=0; ilayer<layerHits_.size(); ++ilayer) {
    dprint("Hits in layer=" << ilayer);
    
#ifdef ETASEG
    segmentMap_[ilayer].resize(Config::nEtaPart);    
    // eta first then phi
    std::sort(layerHits_[ilayer].begin(), layerHits_[ilayer].end(), sortByZ);
    std::vector<int> lay_eta_bin_count(Config::nEtaPart);
    for (int ihit=0;ihit<layerHits_[ilayer].size();++ihit) {
      int etabin = getEtaPartition(layerHits_[ilayer][ihit].eta());
      dprint("ihit: " << ihit << " eta: " << layerHits_[ilayer][ihit].eta() << " etabin: " << etabin);
      lay_eta_bin_count[etabin]++;
    }
    //now set index and size in partitioning map and then sort the bin by phi
    
    int lastEtaIdxFound = -1;
    int lastPhiIdxFound = -1;

    for (int etabin=0; etabin<Config::nEtaPart; ++etabin) {
      int firstEtaBinIdx = lastEtaIdxFound+1;
      int etaBinSize = lay_eta_bin_count[etabin];
      if (etaBinSize>0){
        lastEtaIdxFound+=etaBinSize;
      }

      //sort by phi in each "eta bin"
      std::sort(layerHits_[ilayer].begin() + firstEtaBinIdx,layerHits_[ilayer].begin() + (etaBinSize+firstEtaBinIdx), sortByPhi); // sort from first to last in eta
      std::vector<int> lay_eta_phi_bin_count(Config::nPhiPart);

      for(int ihit = firstEtaBinIdx; ihit < etaBinSize+firstEtaBinIdx; ++ihit){
        dprint("ihit: " << ihit << " r(layer): " << layerHits_[ilayer][ihit].r() << "(" << ilayer << ") phi: " 
	                << layerHits_[ilayer][ihit].phi() << " phipart: " << getPhiPartition(layerHits_[ilayer][ihit].phi()) << " eta: "
	                << layerHits_[ilayer][ihit].eta() << " etapart: " << getEtaPartition(layerHits_[ilayer][ihit].eta()));
        int phibin = getPhiPartition(layerHits_[ilayer][ihit].phi());
        lay_eta_phi_bin_count[phibin]++;
      }

      for (int phibin=0; phibin<Config::nPhiPart; ++phibin) {
        int firstPhiBinIdx = lastPhiIdxFound+1;
        int phiBinSize = lay_eta_phi_bin_count[phibin];
        BinInfo phiBinInfo(firstPhiBinIdx,phiBinSize);
        segmentMap_[ilayer][etabin].push_back(phiBinInfo);
        if (phiBinSize>0){
          lastPhiIdxFound+=phiBinSize;
        }
#ifdef DEBUG
        if ((debug) && (phiBinSize !=0)) dprintf("ilayer: %1u etabin: %1u phibin: %2u first: %2u last: %2u \n", 
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
    std::vector<int> lay_phi_bin_count(Config::nPhiPart);//should it be 63? - yes!
    for (int ihit=0;ihit<layerHits_[ilayer].size();++ihit) {
      dprint("hit r/phi/eta : " << layerHits_[ilayer][ihit].r() << " "
                                << layerHits_[ilayer][ihit].phi() << " " << layerHits_[ilayer][ihit].eta());

      int phibin = getPhiPartition(layerHits_[ilayer][ihit].phi());
      lay_phi_bin_count[phibin]++;
    }

    //now set index and size in partitioning map
    int lastIdxFound = -1;
    for (int bin=0; bin<Config::nPhiPart; ++bin) {
      int binSize = lay_phi_bin_count[bin];
      int firstBinIdx = lastIdxFound+1;
      BinInfo binInfo(firstBinIdx, binSize);
      segmentMap_[ilayer][0].push_back(binInfo); // [0] bin is just the only eta bin ... reduce ifdefs
      if (binSize>0){
        lastIdxFound+=binSize;
      }
    }
#endif
  } // end loop over layers

#ifdef DEBUG
  for (int ilayer = 0; ilayer < Config::nLayers; ilayer++) {
    dmutex_guard;
    int etahitstotal = 0;
    for (int etabin = 0; etabin < Config::nEtaPart; etabin++){
      int etahits = segmentMap_[ilayer][etabin][Config::nPhiPart-1].first + segmentMap_[ilayer][etabin][Config::nPhiPart-1].second - segmentMap_[ilayer][etabin][0].first;
      std::cout << "etabin: " << etabin << " hits in bin: " << etahits << std::endl;
      etahitstotal += etahits;

      for (int phibin = 0; phibin < Config::nPhiPart; phibin++){
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
  buildSeedsByRoadSearch(seedTracks_,seedTracksExtra_,layerHits_,segmentMap_,*this);
  //  buildSeedsByRoadTriplets(seedTracks_,seedTracksExtra_,layerHits_,segmentMap_,*this);
  //buildSeedsByRZFirstRPhiSecond(seedTracks_,seedTracksExtra_,layerHits_,segmentMap_,*this);
#else
  buildSeedsByMC(simTracks_,seedTracks_,seedTracksExtra_,*this);
  simTracksExtra_ = seedTracksExtra_;
#endif
  std::sort(seedTracks_.begin(), seedTracks_.end(), tracksByPhi);
  validation_.alignTrackExtra(seedTracks_,seedTracksExtra_);   // if we sort here, also have to sort seedTracksExtra and redo labels.
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

void Event::Validate(){
  // KM: Config tree just filled once... in main.cc
  if (Config::normal_val) {
    validation_.setTrackExtras(*this);
    validation_.makeSimTkToRecoTksMaps(*this);
    validation_.makeSeedTkToRecoTkMaps(*this);
    validation_.fillEfficiencyTree(*this);
    validation_.fillFakeRateTree(*this);
    if (Config::full_val) {
      validation_.fillSegmentTree(segmentMap_,evtID_);
      validation_.fillBranchTree(evtID_);
      validation_.fillGeometryTree(*this);
      validation_.fillConformalTree(*this);
    }
  }

  if (Config::super_debug) { // super debug mode
    validation_.fillDebugTree(*this);
  }

  if (Config::fit_val) { // fit val for z-phi tuning
    validation_.fillFitTree(*this);
  }
}

void Event::PrintStats(const TrackVec& trks, TrackExtraVec& trkextras)
{
  int miss(0), found(0), fp_10(0), fp_20(0), hit8(0), h8_10(0), h8_20(0);

  for (auto&& trk : trks) {
    auto&& extra = trkextras[trk.label()];
    extra.setMCTrackIDInfo(trk, layerHits_, simHitsInfo_);
    if (extra.mcTrackID() < 0) {
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
 
void Event::write_out(FILE *fp) 
{

  int nt = simTracks_.size();
  fwrite(&nt, sizeof(int), 1, fp);
  fwrite(&simTracks_[0], sizeof(Track), nt, fp);

  if (Config::normal_val || Config::fit_val) {
    for (int it = 0; it<nt; ++it) {
      int nts = simTrackStates_[it].size();
      fwrite(&nts, sizeof(int), 1, fp);
      fwrite(&simTrackStates_[it][0], sizeof(TrackState), nts, fp);
    }
  }

  int nl = layerHits_.size();
  fwrite(&nl, sizeof(int), 1, fp);
  for (int il = 0; il<nl; ++il) {
    int nh = layerHits_[il].size();
    fwrite(&nh, sizeof(int), 1, fp);
    fwrite(&layerHits_[il][0], sizeof(Hit), nh, fp);
  }

  int nm = simHitsInfo_.size();
  fwrite(&nm, sizeof(int), 1, fp);
  fwrite(&simHitsInfo_[0], sizeof(MCHitInfo), nm, fp);

  //layerHitMap_ is recreated afterwards

  /*
  printf("write %i tracks\n",nt);
  for (int it = 0; it<nt; it++) {
    printf("track with pT=%5.3f\n",simTracks_[it].pT());
    for (int ih=0; ih<simTracks_[it].nTotalHits(); ++ih) {
      printf("hit idx=%i\n", simTracks_[it].getHitIdx(ih));
    }
  }
  printf("write %i layers\n",nl);
  for (int il = 0; il<nl; il++) {
    printf("write %i hits in layer %i\n",layerHits_[il].size(),il);
    for (int ih = 0; ih<layerHits_[il].size(); ih++) {
      printf("hit with r=%5.3f x=%5.3f y=%5.3f z=%5.3f\n",layerHits_[il][ih].r(),layerHits_[il][ih].x(),layerHits_[il][ih].y(),layerHits_[il][ih].z());
    }
  }
  */
}

void Event::read_in(FILE *fp)
{

  int nt;
  fread(&nt, sizeof(int), 1, fp);
  simTracks_.resize(nt);
  fread(&simTracks_[0], sizeof(Track), nt, fp);
  Config::nTracks = nt;

  if (Config::normal_val || Config::fit_val) {
    simTrackStates_.resize(nt);
    for (int it = 0; it<nt; ++it) {
      int nts;
      fread(&nts, sizeof(int), 1, fp);
      simTrackStates_[it].resize(nts);
      fread(&simTrackStates_[it][0], sizeof(TrackState), nts, fp);
    }
    // now do the validation copying... ugh
    for (int imc = 0; imc < nt; ++imc) {
      validation_.collectSimTkTSVecMapInfo(imc,simTrackStates_[imc]); 
    }
    simTrackStates_.clear();
  }

  int nl;
  fread(&nl, sizeof(int), 1, fp);
  layerHits_.resize(nl);
  for (int il = 0; il<nl; ++il) {
    int nh;
    fread(&nh, sizeof(int), 1, fp);
    layerHits_[il].resize(nh);
    fread(&layerHits_[il][0], sizeof(Hit), nh, fp);
  }

  int nm; 
  fread(&nm, sizeof(int), 1, fp);
  simHitsInfo_.resize(nm);
  fread(&simHitsInfo_[0], sizeof(MCHitInfo), nm, fp);

  if (Config::useCMSGeom || Config::readCmsswSeeds) {
    int ns;
    fread(&ns, sizeof(int), 1, fp);
    seedTracks_.resize(ns);
    if (Config::readCmsswSeeds) fread(&seedTracks_[0], sizeof(Track), ns, fp);
    else fseek(fp, sizeof(Track)*ns, SEEK_CUR);
    /*
    printf("read %i seedtracks\n",nt);
    for (int it = 0; it<ns; it++) {
      printf("seedtrack with q=%i pT=%5.3f nHits=%i and label=%i\n",seedTracks_[it].charge(),seedTracks_[it].pT(),seedTracks_[it].nFoundHits(),seedTracks_[it].label());
    }
    */
  }

  /*
  printf("read %i simtracks\n",nt);
  for (int it = 0; it<nt; it++) {
    printf("simtrack with q=%i pT=%5.3f and nHits=%i\n",simTracks_[it].charge(),simTracks_[it].pT(),simTracks_[it].nFoundHits());
    for (int ih=0; ih<simTracks_[it].nTotalHits(); ++ih) {
      if (simTracks_[it].getHitIdx(ih)>=0)
	printf("hit #%i idx=%i pos r=%5.3f\n",ih,simTracks_[it].getHitIdx(ih),layerHits_[ih][simTracks_[it].getHitIdx(ih)].r());
      else
	printf("hit #%i idx=%i\n",ih,simTracks_[it].getHitIdx(ih));
    }
  }

  printf("read %i layers\n",nl);
  for (int il = 0; il<nl; il++) {
    printf("read %i hits in layer %i\n",layerHits_[il].size(),il);
    for (int ih = 0; ih<layerHits_[il].size(); ih++) {
      printf("hit with mcHitID=%i r=%5.3f x=%5.3f y=%5.3f z=%5.3f\n",layerHits_[il][ih].mcHitID(),layerHits_[il][ih].r(),layerHits_[il][ih].x(),layerHits_[il][ih].y(),layerHits_[il][ih].z());
    }
  }
  printf("read event done\n",nl);
  */
}
