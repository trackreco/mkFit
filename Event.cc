#include "Event.h"

#include "Simulation.h"
#include "KalmanUtils.h"
#include "seedtest.h"
#include "buildtest.h"
#include "fittest.h"
#include "ConformalUtils.h"

//#define DEBUG
#include "Debug.h"

#ifdef TBB
#include "tbb/tbb.h"
#endif

std::mutex Event::printmutex;

inline bool sortByPhi(const Hit& hit1, const Hit& hit2)
{
  return hit1.phi()<hit2.phi();
}

static bool tracksByPhi(const Track& t1, const Track& t2)
{
  return t1.posPhi()<t2.posPhi();
}

inline bool sortByEta(const Hit& hit1, const Hit& hit2){
  return hit1.eta()<hit2.eta();
}

// within a layer with a "reasonable" geometry, ordering by Z is the same as eta
inline bool sortByZ(const Hit& hit1, const Hit& hit2){
  return hit1.z()<hit2.z();
}

Event::Event(const Geometry& g, Validation& v, int evtID, int threads) :
  geom_(g), validation_(v),
  evtID_(evtID), threads_(threads), mcHitIDCounter_(0)
{
  layerHits_.resize(Config::nTotalLayers);

  validation_.resetValidationMaps(); // need to reset maps for every event.
}

void Event::Reset(int evtID)
{
  evtID_ = evtID;
  mcHitIDCounter_ = 0;

  for (auto&& l : layerHits_) { l.clear(); }

  simHitsInfo_.clear();
  simTrackStates_.clear();
  simTracks_.clear();
  simTracksExtra_.clear();
  seedTracks_.clear();
  seedTracksExtra_.clear();
  candidateTracks_.clear();
  candidateTracksExtra_.clear();
  fitTracks_.clear();
  fitTracksExtra_.clear();

  validation_.resetValidationMaps(); // need to reset maps for every event.
}

void Event::RemapHits(TrackVec & tracks)
{
  std::unordered_map<int,int> simHitMap;
  int max_layer = Config::nTotalLayers;
  for (int ilayer = 0; ilayer < max_layer; ++ilayer)
  {
    const auto & hit_vec = layerHits_[ilayer];
    const auto   size = hit_vec.size();
    for (int index = 0; index < size; ++index)
    {
      simHitMap[hit_vec[index].mcHitID()] = index;
    }
  }
  for (auto&& track : tracks)
  {
    for (int i = 0; i < track.nTotalHits(); ++i)
    {
      int hitidx = track.getHitIdx(i);
      int hitlyr = track.getHitLyr(i);
      if (hitidx >= 0)
      {
        track.setHitIdx(i, simHitMap[layerHits_[hitlyr][hitidx].mcHitID()]);
      }
    }
  }
}

void Event::Simulate()
{
  simTracks_.resize(Config::nTracks);
  simHitsInfo_.reserve(Config::nTotHit * Config::nTracks);
  simTrackStates_.reserve(Config::nTotHit * Config::nTracks);

  for (auto&& l : layerHits_) {
    l.clear();
    l.reserve(Config::nTracks);
  }

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
      MCHitInfoVec hitinfos;
      TSVec  initialTSs;
      // int starting_layer  = 0; --> for displaced tracks, may want to consider running a separate Simulate() block with extra parameters

      int q=0;//set it in setup function
      // do the simulation
      if (Config::useCMSGeom) setupTrackFromTextFile(pos,mom,covtrk,hits,*this,itrack,q,tmpgeom,initialTSs);
      else if (Config::endcapTest) setupTrackByToyMCEndcap(pos,mom,covtrk,hits,*this,itrack,q,tmpgeom,initialTSs);
      else setupTrackByToyMC(pos,mom,covtrk,hits,hitinfos,*this,itrack,q,tmpgeom,initialTSs); 

      // XXMT4K What genius set up separate setupTrackByToyMC / setupTrackByToyMCEndcap?
      // setupTrackByToyMCEndcap does not have scattering and who knows what else.
      // In the original commit Giuseppe he only did fittest ... strangle ... etc ...
      // I'll just review/fix setupTrackByToyMC() for now (so Config::endcapTest = false).
      // See also ifdef just below:

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
      // MT: I'm putting in a mutex for now ...
      std::lock_guard<std::mutex> lock(mcGatherMutex_);

      simTracks_[itrack] = Track(q,pos,mom,covtrk,0.0f);
      auto& sim_track = simTracks_[itrack];
      sim_track.setLabel(itrack);

      // XXKM4MT
      // Sorta assumes one hit per layer -- could just make this 2X nMaxSimHits inside track object (without loopers)
      // This really would only matter for validation and seeding...
      // Could imagine making an inherited class for sim tracks that keeps tracks overlaps
      assert(hits.size() == hitinfos.size());
      for (int i = 0; i < hits.size(); ++i)
      {
        // set to the correct hit index after sorting
        sim_track.addHitIdx(layerHits_[hitinfos[i].layer_].size(), hitinfos[i].layer_, 0.0f);
        layerHits_[hitinfos[i].layer_].emplace_back(hits[i]);

        simHitsInfo_.emplace_back(hitinfos[i]);
	if (Config::root_val || Config::fit_val) 
        {
	  simTrackStates_.emplace_back(initialTSs[i]);
	}
      }
    }
#ifdef TBB
  });
#endif

  // do some work to ensure everything is aligned after multithreading 	
  std::unordered_map<int,int> mcHitIDMap;
  for (int ihit = 0; ihit < simHitsInfo_.size(); ihit++)
  {
    mcHitIDMap[simHitsInfo_[ihit].mcHitID()] = ihit;
  }

  std::sort(simHitsInfo_.begin(),   simHitsInfo_.end(),
	    [](const MCHitInfo& a, const MCHitInfo& b)
	    { return a.mcHitID() < b.mcHitID(); });
		
  TSVec tmpTSVec(simTrackStates_.size());
  for (int its = 0; its < simTrackStates_.size(); its++)
  {
    tmpTSVec[its] = simTrackStates_[mcHitIDMap[its]];
  }
  simTrackStates_ = tmpTSVec;
}

void Event::Segment(BinInfoMap & segmentMap)
{
#ifdef DEBUG
  bool debug=true;
#endif
  segmentMap.resize(Config::nTotalLayers);

  //sort in phi and dump hits per layer, fill phi partitioning
  for (int ilayer=0; ilayer<layerHits_.size(); ++ilayer) {
    dprint("Hits in layer=" << ilayer);
    
    segmentMap[ilayer].resize(Config::nEtaPart);    
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
        segmentMap[ilayer][etabin].push_back(phiBinInfo);
        if (phiBinSize>0){
          lastPhiIdxFound+=phiBinSize;
        }
#ifdef DEBUG
        if ((debug) && (phiBinSize !=0)) dprintf("ilayer: %1u etabin: %1u phibin: %2u first: %2u last: %2u \n", 
                                                 ilayer, etabin, phibin, 
                                                 segmentMap[ilayer][etabin][phibin].first, 
                                                 segmentMap[ilayer][etabin][phibin].second+segmentMap[ilayer][etabin][phibin].first
                                                 );
#endif
      } // end loop over storing phi index
    } // end loop over storing eta index
  } // end loop over layers

#ifdef DEBUG
  for (int ilayer = 0; ilayer < Config::nLayers; ilayer++) {
    dmutex_guard;
    int etahitstotal = 0;
    for (int etabin = 0; etabin < Config::nEtaPart; etabin++){
      int etahits = segmentMap[ilayer][etabin][Config::nPhiPart-1].first + segmentMap[ilayer][etabin][Config::nPhiPart-1].second - segmentMap[ilayer][etabin][0].first;
      std::cout << "etabin: " << etabin << " hits in bin: " << etahits << std::endl;
      etahitstotal += etahits;

      for (int phibin = 0; phibin < Config::nPhiPart; phibin++){
	//	if (segmentMap[ilayer][etabin][phibin].second > 3) {std::cout << "   phibin: " << phibin << " hits: " << segmentMap[ilayer][etabin][phibin].second << std::endl;}
      }
    }
    std::cout << "layer: " << ilayer << " totalhits: " << etahitstotal << std::endl;
  }
#endif

  // need to reset simtrack hit indices after sorting!
  RemapHits(simTracks_);
}

void Event::Seed(const BinInfoMap & segmentMap)
{
#ifdef ENDTOEND
  buildSeedsByRoadSearch(seedTracks_,seedTracksExtra_,layerHits_,segmentMap,*this);
  //  buildSeedsByRoadTriplets(seedTracks_,seedTracksExtra_,layerHits_,segmentMap,*this);
  //buildSeedsByRZFirstRPhiSecond(seedTracks_,seedTracksExtra_,layerHits_,segmentMap,*this);
#else
  buildSeedsByMC(simTracks_,seedTracks_,seedTracksExtra_,*this);
  simTracksExtra_ = seedTracksExtra_;
#endif
  std::sort(seedTracks_.begin(), seedTracks_.end(), tracksByPhi);
  validation_.alignTrackExtra(seedTracks_,seedTracksExtra_);   // if we sort here, also have to sort seedTracksExtra and redo labels.
}

void Event::Find(const BinInfoMap & segmentMap)
{
  buildTracksBySeeds(segmentMap,*this);
  //  buildTracksByLayers(segmentMap,*this);

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

void Event::Validate()
{
  // KM: Config tree just filled once... in main.cc
  if (Config::root_val) {
    validation_.setTrackExtras(*this);
    validation_.makeSimTkToRecoTksMaps(*this);
    validation_.makeSeedTkToRecoTkMaps(*this);
    validation_.fillEfficiencyTree(*this);
    validation_.fillFakeRateTree(*this);
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
    extra.setMCTrackIDInfoByLabel(trk, layerHits_, simHitsInfo_);
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
  static std::mutex writemutex;
  std::lock_guard<std::mutex> writelock(writemutex);

  auto start = ftell(fp);
  int evsize = sizeof(int);
  fwrite(&evsize, sizeof(int), 1, fp); // this will be overwritten at the end

  int nt = simTracks_.size();
  fwrite(&nt, sizeof(int), 1, fp);
  fwrite(&simTracks_[0], sizeof(Track), nt, fp);
  evsize += sizeof(int) + nt*sizeof(Track);

  if (Config::root_val || Config::fit_val) {
    int nts = simTrackStates_.size();
    fwrite(&nts, sizeof(int), 1, fp);
    fwrite(&simTrackStates_[0], sizeof(TrackState), nts, fp);
    evsize += sizeof(int) + nts*sizeof(TrackState);
  }

  int nl = layerHits_.size();
  fwrite(&nl, sizeof(int), 1, fp);
  evsize += sizeof(int);
  for (int il = 0; il<nl; ++il) {
    int nh = layerHits_[il].size();
    fwrite(&nh, sizeof(int), 1, fp);
    fwrite(&layerHits_[il][0], sizeof(Hit), nh, fp);
    evsize += sizeof(int) + nh*sizeof(Hit);
  }

  int nm = simHitsInfo_.size();
  fwrite(&nm, sizeof(int), 1, fp);
  fwrite(&simHitsInfo_[0], sizeof(MCHitInfo), nm, fp);
  evsize += sizeof(int) + nm*sizeof(MCHitInfo);

  if (Config::useCMSGeom || Config::readCmsswSeeds) {
    int ns = seedTracks_.size();
    fwrite(&ns, sizeof(int), 1, fp);
    fwrite(&seedTracks_[0], sizeof(Track), ns, fp);
    evsize += sizeof(int) + ns*sizeof(Track);
  }

  fseek(fp, start, SEEK_SET);
  fwrite(&evsize, sizeof(int), 1, fp);
  fseek(fp, 0, SEEK_END);
  
  //layerHitMap_ is recreated afterwards

  /*
  printf("write %i tracks\n",nt);
  for (int it = 0; it<nt; it++) {
    printf("track with pT=%5.3f\n",simTracks_[it].pT());
    for (int ih=0; ih<simTracks_[it].nTotalHits(); ++ih) {
      printf("hit lyr:%2d idx=%i\n", simTracks_[it].getHitLyr(ih), simTracks_[it].getHitIdx(ih));
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

void Event::read_in(FILE *fp, int version)
{
  static long pos = sizeof(int); // header size
  int evsize;

  if (version > 0) {
    static std::mutex readmutex;
    std::lock_guard<std::mutex> readlock(readmutex);

    fseek(fp, pos, SEEK_SET);
    fread(&evsize, sizeof(int), 1, fp);
    pos += evsize;
  }

  int nt;
  fread(&nt, sizeof(int), 1, fp);
  simTracks_.resize(nt);
  fread(&simTracks_[0], sizeof(Track), nt, fp);
  Config::nTracks = nt;

  if (Config::root_val || Config::fit_val)
  {
    int nts; 
    fread(&nts, sizeof(int), 1, fp);
    simTrackStates_.resize(nts);
    fread(&simTrackStates_[0], sizeof(TrackState), nts, fp);
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

  /* */
  printf("read %i simtracks\n",nt);
  /*
  for (int it = 0; it<nt; it++) {
    printf("simtrack %d with q=%i pT=%5.3f eta=%6.3f and nHits=%i\n",it,simTracks_[it].charge(),simTracks_[it].pT(),simTracks_[it].momEta(),simTracks_[it].nFoundHits());
    for (int ih=0; ih<simTracks_[it].nTotalHits(); ++ih) {
      int lyr = simTracks_[it].getHitLyr(ih);
      int idx = simTracks_[it].getHitIdx(ih);
      if (idx >= 0)
	printf("  hit #%i lyr=%d idx=%i pos r=%5.3f z=%6.3f\n",
               ih, lyr, idx, layerHits_[lyr][idx].r(), layerHits_[lyr][idx].z());
      else
	printf("  hit #%i idx=%i\n",ih,simTracks_[it].getHitIdx(ih));
    }
  }
  printf("read %i layers\n",nl);
  int total_hits = 0;
  for (int il = 0; il<nl; il++) {
    printf("read %i hits in layer %i\n",layerHits_[il].size(),il);
    total_hits += layerHits_[il].size();
    for (int ih = 0; ih<layerHits_[il].size(); ih++) {
      printf("  hit with mcHitID=%i r=%5.3f x=%5.3f y=%5.3f z=%5.3f\n",layerHits_[il][ih].mcHitID(),layerHits_[il][ih].r(),layerHits_[il][ih].x(),layerHits_[il][ih].y(),layerHits_[il][ih].z());
    }
  }
  printf("total_hits = %d\n", total_hits);
  */
  printf("read event done\n");
}
