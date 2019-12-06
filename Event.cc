#include "Event.h"

#include "Simulation.h"
#include "KalmanUtils.h"
#include "seedtest.h"
#include "buildtest.h"
#include "fittest.h"
#include "ConformalUtils.h"
#include "TrackerInfo.h"

//#define DEBUG
#include "Debug.h"

#ifdef TBB
#include "tbb/tbb.h"
#endif

#include <memory>

namespace
{
  mkfit::Geometry dummyGeometry;
  std::unique_ptr<mkfit::Validation> dummyValidation( mkfit::Validation::make_validation("dummy") );
}

namespace mkfit {

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

void Event::reset_nan_n_silly_counters()
{
  nan_n_silly_per_layer_count_ = 0;
}

Event::Event(int evtID) :
  geom_(dummyGeometry), validation_(*dummyValidation),
  evtID_(evtID), threads_(1), mcHitIDCounter_(0)
{
  reset_nan_n_silly_counters();
  layerHits_.resize(Config::nTotalLayers);
}

Event::Event(const Geometry& g, Validation& v, int evtID, int threads) :
  geom_(g), validation_(v),
  evtID_(evtID), threads_(threads), mcHitIDCounter_(0)
{
  reset_nan_n_silly_counters();
  layerHits_.resize(Config::nTotalLayers);

  validation_.resetValidationMaps(); // need to reset maps for every event.
}

void Event::Reset(int evtID)
{
  evtID_ = evtID;
  mcHitIDCounter_ = 0;
  reset_nan_n_silly_counters();

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
  cmsswTracks_.clear();
  cmsswTracksExtra_.clear();

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
    for (size_t index = 0; index < size; ++index)
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
  simHitsInfo_.reserve(Config::nAvgSimHits * Config::nTracks);
  simTrackStates_.reserve(Config::nAvgSimHits * Config::nTracks);

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
      setupTrackByToyMC(pos,mom,covtrk,hits,hitinfos,*this,itrack,q,tmpgeom,initialTSs); 

      // convert from global cartesian to CCS
      float pt = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
      mom=SVector3(1./pt,atan2(mom[1],mom[0]),atan2(pt,mom[2]));
      for (size_t its = 0; its < initialTSs.size(); its++){
	initialTSs[its].convertFromCartesianToCCS();
      }

      // MT: I'm putting in a mutex for now ...
      std::lock_guard<std::mutex> lock(mcGatherMutex_);

      simTracks_[itrack] = Track(q,pos,mom,covtrk,0.0f);
      auto& sim_track = simTracks_[itrack];
      sim_track.setLabel(itrack);

      // XXKM4MT
      // Sorta assumes one hit per layer.
      // This really would only matter for validation and seeding...
      // Could imagine making an inherited class for sim tracks that keeps tracks overlaps
      assert(hits.size() == hitinfos.size());
      for (size_t i = 0; i < hits.size(); ++i)
      {
        // set to the correct hit index after sorting
        sim_track.addHitIdx(layerHits_[hitinfos[i].layer_].size(), hitinfos[i].layer_, 0.0f);
        layerHits_[hitinfos[i].layer_].emplace_back(hits[i]);

        simHitsInfo_.emplace_back(hitinfos[i]);
	if (Config::sim_val || Config::fit_val) 
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
  for (size_t ihit = 0; ihit < simHitsInfo_.size(); ihit++)
  {
    mcHitIDMap[simHitsInfo_[ihit].mcHitID()] = ihit;
  }

  std::sort(simHitsInfo_.begin(),   simHitsInfo_.end(),
	    [](const MCHitInfo& a, const MCHitInfo& b)
	    { return a.mcHitID() < b.mcHitID(); });
		
  TSVec tmpTSVec(simTrackStates_.size());
  for (size_t its = 0; its < simTrackStates_.size(); its++)
  {
    tmpTSVec[its] = simTrackStates_[mcHitIDMap[its]];
  }
  simTrackStates_.swap(tmpTSVec);
}

void Event::Segment(BinInfoMap & segmentMap)
{
#ifdef DEBUG
  bool debug=true;
#endif
  segmentMap.resize(Config::nTotalLayers);

  //sort in phi and dump hits per layer, fill phi partitioning
  for (size_t ilayer=0; ilayer<layerHits_.size(); ++ilayer) {
    dprint("Hits in layer=" << ilayer);
    
    segmentMap[ilayer].resize(Config::nEtaPart);    
    // eta first then phi
    std::sort(layerHits_[ilayer].begin(), layerHits_[ilayer].end(), sortByZ);
    std::vector<int> lay_eta_bin_count(Config::nEtaPart);
    for (size_t ihit=0;ihit<layerHits_[ilayer].size();++ihit) {
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
  validation_.alignTracks(seedTracks_,seedTracksExtra_,true);   // if we sort here, also have to sort seedTracksExtra and redo labels.
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
  // special map needed for sim_val_for_cmssw + set the track scores
  if (Config::sim_val_for_cmssw) {
    validation_.makeRecoTkToSeedTkMapsDumbCMSSW(*this);
    validation_.setTrackScoresDumbCMSSW(*this);
  }

  // standard eff/fr/dr validation
  if (Config::sim_val || Config::sim_val_for_cmssw) {
    validation_.setTrackExtras(*this);
    validation_.makeSimTkToRecoTksMaps(*this);
    validation_.makeSeedTkToRecoTkMaps(*this);
    validation_.fillEfficiencyTree(*this);
    validation_.fillFakeRateTree(*this);
  }

  // special cmssw to mkfit validation
  if (Config::cmssw_val) {
    validation_.makeCMSSWTkToSeedTkMap(*this);
    validation_.makeRecoTkToRecoTkMaps(*this);
    validation_.setTrackExtras(*this);
    validation_.makeCMSSWTkToRecoTksMaps(*this);
    validation_.fillCMSSWEfficiencyTree(*this);
    validation_.fillCMSSWFakeRateTree(*this);
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
    extra.setMCTrackIDInfo(trk, layerHits_, simHitsInfo_, simTracks_, false, true);
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
 
void Event::write_out(DataFile &data_file)
{
  FILE *fp = data_file.f_fp;

  static std::mutex writemutex;
  std::lock_guard<std::mutex> writelock(writemutex);

  auto start = ftell(fp);
  int evsize = sizeof(int);
  fwrite(&evsize, sizeof(int), 1, fp); // this will be overwritten at the end

  evsize += write_tracks(fp, simTracks_);

  if (data_file.HasSimTrackStates()) {
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

  if (data_file.HasSeeds()) {
    evsize += write_tracks(fp, seedTracks_);
  }

  if (data_file.HasCmsswTracks()) {
    evsize += write_tracks(fp, cmsswTracks_);
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

// #define DUMP_SEEDS
// #define DUMP_SEED_HITS
// #define DUMP_TRACKS
// #define DUMP_TRACK_HITS
// #define DUMP_LAYER_HITS
// #define DUMP_REC_TRACKS
// #define DUMP_REC_TRACK_HITS

void Event::read_in(DataFile &data_file, FILE *in_fp)
{
  FILE *fp = in_fp ? in_fp : data_file.f_fp;

  data_file.AdvancePosToNextEvent(fp);

  int nt = read_tracks(fp, simTracks_);
  Config::nTracks = nt;

  if (data_file.HasSimTrackStates())
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

  if (data_file.HasSeeds()) {
    int ns = read_tracks(fp, seedTracks_, Config::seedInput != cmsswSeeds);
    (void) ns;

#ifdef DUMP_SEEDS
    printf("Read %i seedtracks (neg value means actual reading was skipped)\n", ns);
    for (int it = 0; it < ns; it++)
    {
      const Track& ss = seedTracks_[it];
      printf("  %3i q=%+i pT=%7.3f eta=% 7.3f nHits=%i label=% i\n",
             it,ss.charge(),ss.pT(),ss.momEta(),ss.nFoundHits(),ss.label());
#ifdef DUMP_SEED_HITS
      for (int ih = 0; ih < seedTracks_[it].nTotalHits(); ++ih)
      {
        int lyr = seedTracks_[it].getHitLyr(ih);
        int idx = seedTracks_[it].getHitIdx(ih);
        if (idx >= 0)
        {
          const Hit &hit = layerHits_[lyr][idx];
          printf("    hit %2d lyr=%3d idx=%4d pos r=%7.3f z=% 8.3f   mc_hit=%3d mc_trk=%3d\n",
                 ih, lyr, idx, layerHits_[lyr][idx].r(), layerHits_[lyr][idx].z(),
                 hit.mcHitID(), hit.mcTrackID(simHitsInfo_));
        }
        else
          printf("    hit %2d idx=%i\n",ih,seedTracks_[it].getHitIdx(ih));

      }
#endif
    }
#endif
  }

  int nert = -99999;
  if (data_file.HasCmsswTracks())
  {
    nert = read_tracks(fp, cmsswTracks_, ! Config::readCmsswTracks);
    (void) nert;
  }

  /*
    // HACK TO ONLY SELECT ONE PROBLEMATIC TRACK.
    // Note that MC matching gets screwed.
    // Works for MC seeding.
    //
    printf("************** SIM SELECTION HACK IN FORCE ********************\n");
    TrackVec x;
    x.push_back(simTracks_[3]);
    simTracks_.swap(x);
    nt = 1;
  */

#ifdef DUMP_TRACKS
  printf("Read %i simtracks\n", nt);
  for (int it = 0; it < nt; it++)
  {
    const Track &t = simTracks_[it];
    printf("  %3i q=%+i pT=%7.3f eta=% 7.3f nHits=%2d  label=%4d\n",
           it, t.charge(), t.pT(), t.momEta(), t.nFoundHits(), t.label());
#ifdef DUMP_TRACK_HITS
    for (int ih = 0; ih < t.nTotalHits(); ++ih)
    {
      int lyr = t.getHitLyr(ih);
      int idx = t.getHitIdx(ih);
      if (idx >= 0)
      {
        const Hit &hit = layerHits_[lyr][idx];
	printf("    hit %2d lyr=%2d idx=%3d pos r=%7.3f x=% 8.3f y=% 8.3f z=% 8.3f   mc_hit=%3d mc_trk=%3d\n",
               ih, lyr, idx, layerHits_[lyr][idx].r(), layerHits_[lyr][idx].x(), layerHits_[lyr][idx].y(), layerHits_[lyr][idx].z(),
               hit.mcHitID(), hit.mcTrackID(simHitsInfo_));
      }
      else
	printf("    hit %2d idx=%i\n", ih, t.getHitIdx(ih));
    }
#endif
  }
#endif
#ifdef DUMP_LAYER_HITS
  printf("Read %i layers\n",nl);
  int total_hits = 0;
  for (int il = 0; il < nl; il++)
  {
    if (layerHits_[il].empty()) continue;

    printf("Read %i hits in layer %i\n", (int) layerHits_[il].size(), il);
    total_hits += layerHits_[il].size();
    for (int ih = 0; ih < (int) layerHits_[il].size(); ih++)
    {
      const Hit &hit = layerHits_[il][ih];
      printf("  mcHitID=%5d r=%10g x=%10g y=%10g z=%10g  sx=%10.4g sy=%10.4e sz=%10.4e\n",
             hit.mcHitID(), hit.r(), hit.x(), hit.y(), hit.z(),
             std::sqrt(hit.exx()), std::sqrt(hit.eyy()), std::sqrt(hit.ezz()));
    }
  }
  printf("Total hits in all layers = %d\n", total_hits);
#endif
#ifdef DUMP_REC_TRACKS
  printf("Read %i rectracks\n", nert);
  for (int it = 0; it < nert; it++)
  {
    const Track &t = cmsswTracks_[it];
    printf("  %i with q=%+i pT=%7.3f eta=% 7.3f nHits=%2d  label=%4d\n",
           it, t.charge(), t.pT(), t.momEta(), t.nFoundHits(), t.label());
#ifdef DUMP_REC_TRACK_HITS
    for (int ih = 0; ih < t.nTotalHits(); ++ih)
    {
      int lyr = t.getHitLyr(ih);
      int idx = t.getHitIdx(ih);
      if (idx >= 0)
      {
        const Hit &hit = layerHits_[lyr][idx];
        printf("    hit %2d lyr=%2d idx=%3d pos r=%7.3f z=% 8.3f   mc_hit=%3d mc_trk=%3d\n",
               ih, lyr, idx, hit.r(), hit.z(),
               hit.mcHitID(), hit.mcTrackID(simHitsInfo_));
      }
      else
       printf("    hit %2d        idx=%i\n", ih, t.getHitIdx(ih));
    }
#endif
  }
#endif

  if (Config::kludgeCmsHitErrors)
  {
    kludge_cms_hit_errors();
  }

  if (!Config::silent) printf("Read complete, %d simtracks on file.\n", nt);
}

//------------------------------------------------------------------------------

int Event::write_tracks(FILE *fp, const TrackVec& tracks)
{
  // Returns total number of bytes written.

  int n_tracks  = tracks.size();
  fwrite(&n_tracks, sizeof(int), 1, fp);

  auto start = ftell(fp);
  int data_size = 2 * sizeof(int) + n_tracks * sizeof(Track);
  fwrite(&data_size, sizeof(int), 1, fp);

  fwrite(tracks.data(), sizeof(Track), n_tracks, fp);

  for (int i = 0; i < n_tracks; ++i)
  {
    fwrite(tracks[i].BeginHitsOnTrack(), sizeof(HitOnTrack), tracks[i].nTotalHits(), fp);
    data_size += tracks[i].nTotalHits() * sizeof(HitOnTrack);
  }

  fseek(fp, start, SEEK_SET);
  fwrite(&data_size, sizeof(int), 1, fp);
  fseek(fp, 0, SEEK_END);

  return data_size;
}

int Event::read_tracks(FILE *fp, TrackVec& tracks, bool skip_reading)
{
  // Returns number of read tracks (negative if actual reading was skipped).

  int n_tracks, data_size;
  fread(&n_tracks,  sizeof(int), 1, fp);
  fread(&data_size, sizeof(int), 1, fp);

  if (skip_reading)
  {
    fseek(fp, data_size, SEEK_CUR);
    n_tracks = -n_tracks;
  }
  else
  {
    tracks.resize(n_tracks);

    fread(tracks.data(), sizeof(Track), n_tracks, fp);

    for (int i = 0; i < n_tracks; ++i)
    {
      tracks[i].resizeHitsForInput();
      fread(tracks[i].BeginHitsOnTrack_nc(), sizeof(HitOnTrack), tracks[i].nTotalHits(), fp);
    }
  }

  return n_tracks;
}

//------------------------------------------------------------------------------

void Event::setInputFromCMSSW(std::vector<HitVec> hits, TrackVec seeds)
{
  layerHits_ = std::move(hits);
  seedTracks_ = std::move(seeds);
}

//------------------------------------------------------------------------------

void Event::kludge_cms_hit_errors()
{
  // Enforce Vxy on all layers, Vz on pixb only.

  const float Exy = 15 * 1e-4, Vxy = Exy * Exy;
  const float Ez  = 30 * 1e-4, Vz  = Ez  * Ez;

  int nl = layerHits_.size();

  int cnt = 0;

  for (int il = 0; il < nl; il++)
  {
    if (layerHits_[il].empty()) continue;

    for (Hit & h : layerHits_[il])
    {
      SVector6 &c = h.error_nc();

      float vxy = c[0] + c[2];
      if (vxy < Vxy)
      {
        c[0] *= Vxy / vxy;
        c[2] *= Vxy / vxy;
        ++cnt;
      }
      if (il < 4 && c[5] < Vz)
      {
        c[5] = Vz;
        ++cnt;
      }
    }
  }

  printf("Event::kludge_cms_hit_errors processed %d layers, kludged %d entries.\n", nl, cnt);
}

//------------------------------------------------------------------------------

int Event::clean_cms_simtracks()
{
  // Sim tracks from cmssw have the following issues:
  // - hits are not sorted by layer;
  // - there are tracks with too low number of hits, even 0;
  // - even with enough hits, there can be too few layers (esp. in endcap);
  // - tracks from secondaries can have extremely low pT.
  // Possible further checks:
  // - make sure enough hits exist in seeding layers.
  //
  // What is done:
  // 1. Hits are sorted by layer;
  // 2. Non-findable tracks are marked with Track::Status::not_findable flag.
  //
  // Returns number of passed simtracks.

  dprintf("Event::clean_cms_simtracks processing %d simtracks.\n", simTracks_.size());

  int n_acc = 0;
  int i = -1;//wrap in ifdef DEBUG?
  for (Track & t : simTracks_)
  {
    i++;

    t.sortHitsByLayer();

    const int lyr_cnt = t.nUniqueLayers();

    //const int lasthit = t.getLastFoundHitPos();
    //const float eta = layerHits_[t.getHitLyr(lasthit)][t.getHitIdx(lasthit)].eta();

    if (lyr_cnt < Config::cmsSelMinLayers) // || Config::TrkInfo.is_transition(eta))
    {
      dprintf("Rejecting simtrack %d, n_hits=%d, n_layers=%d, pT=%f\n", i, t.nFoundHits(), lyr_cnt, t.pT());
      t.setNotFindable();
    }
    else
    {
      dprintf("Accepting simtrack %d, n_hits=%d, n_layers=%d, pT=%f\n", i, t.nFoundHits(), lyr_cnt, t.pT());
      ++n_acc;
    }
  }

  return n_acc;
}

void Event::print_tracks(const TrackVec& tracks, bool print_hits) const
{
  const int nt = tracks.size();

  printf("Event::print_tracks printing %d tracks %s hits:\n", nt, (print_hits ? "with" : "without"));
  for (int it = 0; it < nt; it++)
  {
    const Track &t = tracks[it];
    printf("  %i with q=%+i pT=%7.3f eta=% 7.3f nHits=%2d  label=%4d findable=%d\n",
           it, t.charge(), t.pT(), t.momEta(), t.nFoundHits(), t.label(), t.isFindable());

    if (print_hits)
    {
      for (int ih = 0; ih < t.nTotalHits(); ++ih)
      {
        int lyr = t.getHitLyr(ih);
        int idx = t.getHitIdx(ih);
        if (idx >= 0)
        {
          const Hit &hit = layerHits_[lyr][idx];
          printf("    hit %2d lyr=%2d idx=%3d pos r=%7.3f z=% 8.3f   mc_hit=%3d mc_trk=%3d\n",
                 ih, lyr, idx, layerHits_[lyr][idx].r(), layerHits_[lyr][idx].z(),
                 hit.mcHitID(), hit.mcTrackID(simHitsInfo_));
        }
        else
          printf("    hit %2d lyr=%2d idx=%3d\n", ih, t.getHitLyr(ih), t.getHitIdx(ih));
      }
    }
  }
}

int Event::clean_cms_seedtracks()
{

  const float etamax_brl = Config::c_etamax_brl;
  const float dpt_brl_0  = Config::c_dpt_brl_0;
  const float dpt_ec_0   = Config::c_dpt_ec_0;
  const float ptmax_0    = Config::c_ptmax_0;
  const float dpt_1      = Config::c_dpt_1;
  const float ptmax_1    = Config::c_ptmax_1;
  const float dpt_2      = Config::c_dpt_2;
  const float ptmax_2    = Config::c_ptmax_2;
  const float dpt_3      = Config::c_dpt_3;
  const float dzmax_brl  = Config::c_dzmax_brl;
  const float drmax_brl  = Config::c_drmax_brl;
  const float ptmin_hpt  = Config::c_ptmin_hpt;
  const float dzmax_hpt  = Config::c_dzmax_hpt;
  const float drmax_hpt  = Config::c_drmax_hpt;
  const float dzmax_els  = Config::c_dzmax_els;
  const float drmax_els  = Config::c_drmax_els;

  const float dzmax2_brl = dzmax_brl*dzmax_brl;
  const float drmax2_brl = drmax_brl*drmax_brl;
  const float dzmax2_hpt = dzmax_hpt*dzmax_hpt;
  const float drmax2_hpt = drmax_hpt*drmax_hpt;
  const float dzmax2_els = dzmax_els*dzmax_els;
  const float drmax2_els = drmax_els*drmax_els;

  const int ns = seedTracks_.size();

  TrackVec cleanSeedTracks;
  cleanSeedTracks.reserve(ns);
  std::vector<bool> writetrack(ns, true);

  const float invR1GeV = 1.f/Config::track1GeVradius;

  std::vector<int>    nHits(ns);
  std::vector<int>    charge(ns);
  std::vector<float>  oldPhi(ns);
  std::vector<float>  pos2(ns);
  std::vector<float>  eta(ns);
  std::vector<float>  theta(ns);
  std::vector<float>  invptq(ns);
  std::vector<float>  pt(ns);
  std::vector<float>  x(ns);
  std::vector<float>  y(ns);
  std::vector<float>  z(ns);

  for(int ts=0; ts<ns; ts++){
    const Track & tk = seedTracks_[ts];
    nHits[ts] = tk.nFoundHits();
    charge[ts] = tk.charge();
    oldPhi[ts] = tk.momPhi();
    pos2[ts] = std::pow(tk.x(), 2) + std::pow(tk.y(), 2);
    eta[ts] = tk.momEta();
    theta[ts] = std::atan2(tk.pT(),tk.pz());
    invptq[ts] = tk.charge()*tk.invpT();
    pt[ts] = tk.pT();
    x[ts] = tk.x();
    y[ts] = tk.y();
    z[ts] = tk.z();
  }

  for(int ts=0; ts<ns; ts++){

    if (not writetrack[ts]) continue;//FIXME: this speed up prevents transitive masking; check build cost!

    const float oldPhi1 = oldPhi[ts];
    const float pos2_first = pos2[ts];
    const float Eta1 = eta[ts];
    const float Pt1 = pt[ts];
    const float invptq_first = invptq[ts]; 

    //#pragma simd /* Vectorization via simd had issues with icc */
    for (int tss= ts+1; tss<ns; tss++){

      const float Pt2 = pt[tss];

      ////// Always require charge consistency. If different charge is assigned, do not remove seed-track
      if(charge[tss] != charge[ts])
        continue;
      
      const float thisDPt = std::abs(Pt2-Pt1);
      ////// Require pT consistency between seeds. If dpT is large, do not remove seed-track.
      ////// Adaptive thresholds, based on pT of reference seed-track (choice is a compromise between efficiency and duplicate rate):
      ////// - 2.5% if track is barrel and w/ pT<2 GeV
      ////// - 1.25% if track is non-barrel and w/ pT<2 GeV
      ////// - 10% if track w/ 2<pT<5 GeV
      ////// - 20% if track w/ 5<pT<10 GeV
      ////// - 25% if track w/ pT>10 GeV
      if(thisDPt>dpt_brl_0*(Pt1) && Pt1<ptmax_0 && std::abs(Eta1)<etamax_brl)
	continue;

      else if(thisDPt>dpt_ec_0*(Pt1) && Pt1<ptmax_0 && std::abs(Eta1)>etamax_brl)
	continue;

      else if(thisDPt>dpt_1*(Pt1) && Pt1>ptmax_0 && Pt1<ptmax_1)
	continue;

      else if(thisDPt>dpt_2*(Pt1) && Pt1>ptmax_1 && Pt1<ptmax_2)
	continue;

      else if(thisDPt>dpt_3*(Pt1) && Pt1>ptmax_2)
	continue;


      const float Eta2 = eta[tss];
      const float deta2 = std::pow(Eta1-Eta2, 2);

      const float oldPhi2 = oldPhi[tss];

      const float pos2_second = pos2[tss];
      const float thisDXYSign05 = pos2_second > pos2_first ? -0.5f : 0.5f;

      const float thisDXY = thisDXYSign05*sqrt( std::pow(x[ts]-x[tss], 2) + std::pow(y[ts]-y[tss], 2) );
      
      const float invptq_second = invptq[tss];

      const float newPhi1 = oldPhi1-thisDXY*invR1GeV*invptq_first;
      const float newPhi2 = oldPhi2+thisDXY*invR1GeV*invptq_second;

      const float dphi = cdist(std::abs(newPhi1-newPhi2));

      const float dr2 = deta2+dphi*dphi;
      
      const float thisDZ = z[ts]-z[tss]-thisDXY*(1.f/std::tan(theta[ts])+1.f/std::tan(theta[tss]));
      const float dz2 = thisDZ*thisDZ;

      ////// Reject tracks within dR-dz elliptical window.
      ////// Adaptive thresholds, based on observation that duplicates are more abundant at large pseudo-rapidity and low track pT
      if(std::abs(Eta1)<etamax_brl){
	if(dz2/dzmax2_brl+dr2/drmax2_brl<1.0f)
	  writetrack[tss]=false;	
      }
      else if(Pt1>ptmin_hpt){
	if(dz2/dzmax2_hpt+dr2/drmax2_hpt<1.0f)
	  writetrack[tss]=false;
      }
      else {
	if(dz2/dzmax2_els+dr2/drmax2_els<1.0f)
	  writetrack[tss]=false;
      }

    }
   
    if(writetrack[ts])
      cleanSeedTracks.emplace_back(seedTracks_[ts]);

  }
  
  seedTracks_.swap(cleanSeedTracks);

#ifdef DEBUG
  {
    const int ns2 = seedTracks_.size();
    printf("Number of CMS seeds before %d --> after %d cleaning\n", ns, ns2);

    for (int it = 0; it < ns2; it++)
    {
      const Track& ss = seedTracks_[it];
      printf("  %3i q=%+i pT=%7.3f eta=% 7.3f nHits=%i label=% i\n",
             it,ss.charge(),ss.pT(),ss.momEta(),ss.nFoundHits(),ss.label());
    }
  }
#endif

  return seedTracks_.size();
}

int Event::clean_cms_seedtracks_badlabel()
{
  printf("***\n*** REMOVING SEEDS WITH BAD LABEL. This is a development hack. ***\n***\n");
  TrackVec buf; seedTracks_.swap(buf);
  std::copy_if(buf.begin(), buf.end(), std::back_inserter(seedTracks_), [](const Track& t){ return t.label() >= 0; });
  return seedTracks_.size();
}

int Event::use_seeds_from_cmsswtracks()
{
  int ns = seedTracks_.size();

  TrackVec cleanSeedTracks;
  cleanSeedTracks.reserve(ns);

  for (auto&& cmsswtrack : cmsswTracks_)
  {
    cleanSeedTracks.emplace_back(seedTracks_[cmsswtrack.label()]);
  }

  seedTracks_.swap(cleanSeedTracks);

  return seedTracks_.size();
}

void Event::relabel_bad_seedtracks()
{
  int newlabel = 0;
  for (auto&& track : seedTracks_)
  { 
    if (track.label() < 0) track.setLabel(--newlabel);
  }
}

void Event::relabel_cmsswtracks_from_seeds()
{
  std::map<int,int> cmsswLabelMap;
  for (size_t iseed = 0; iseed < seedTracks_.size(); iseed++)
  {
    for (size_t icmssw = 0; icmssw < cmsswTracks_.size(); icmssw++)
    {
      if (cmsswTracks_[icmssw].label() == static_cast<int>(iseed))
      {
	cmsswLabelMap[icmssw] = seedTracks_[iseed].label();
  	break;
      }
    }
  }
  for (size_t icmssw = 0; icmssw < cmsswTracks_.size(); icmssw++)
  {
    cmsswTracks_[icmssw].setLabel(cmsswLabelMap[icmssw]);
  }
}

//==============================================================================
// DataFile
//==============================================================================

int DataFile::OpenRead(const std::string& fname, bool set_n_layers)
{
  constexpr int min_ver = 4;
  constexpr int max_ver = 4;

  f_fp = fopen(fname.c_str(), "r");
  assert (f_fp != 0 && "Opening of input file failed.");

  fread(&f_header, sizeof(DataFileHeader), 1, f_fp);

  if (f_header.f_magic != 0xBEEF)
  {
    fprintf(stderr, "Incompatible input file (wrong magick).\n");
    exit(1);
  }
  if (f_header.f_format_version < min_ver || f_header.f_format_version > max_ver)
  {
    fprintf(stderr, "Unsupported file version %d. Supported versions are from %d to %d.\n",
            f_header.f_format_version, min_ver, max_ver);
    exit(1);
  }
  if (f_header.f_sizeof_track != sizeof(Track))
  {
    fprintf(stderr, "sizeof(Track) on file (%d) different from current value (%d).\n",
            f_header.f_sizeof_track, (int) sizeof(Track));
    exit(1);
  }
  if (f_header.f_sizeof_hit != sizeof(Hit))
  {
    fprintf(stderr, "sizeof(Hit) on file (%d) different from current value (%d).\n",
            f_header.f_sizeof_hit, (int) sizeof(Hit));
    exit(1);
  }
  if (f_header.f_sizeof_hot != sizeof(HitOnTrack))
  {
    fprintf(stderr, "sizeof(HitOnTrack) on file (%d) different from current value (%d).\n",
            f_header.f_sizeof_hot, (int) sizeof(HitOnTrack));
    exit(1);
  }
  if (set_n_layers)
  {
    Config::nTotalLayers = f_header.f_n_layers;
  }
  else if (f_header.f_n_layers != Config::nTotalLayers)
  {
    fprintf(stderr, "Number of layers on file (%d) is different from current Config::nTotalLayers (%d).\n",
            f_header.f_n_layers, Config::nTotalLayers);
    exit(1);
  }

  printf("Opened file '%s', format version %d, n_layers %d, n_events %d\n",
         fname.c_str(), f_header.f_format_version, f_header.f_n_layers, f_header.f_n_events);
  if (f_header.f_extra_sections)
  {
    printf("  Extra sections:");
    if (f_header.f_extra_sections & ES_SimTrackStates) printf(" SimTrackStates");
    if (f_header.f_extra_sections & ES_Seeds)          printf(" Seeds");
    if (f_header.f_extra_sections & ES_CmsswTracks)    printf(" CmsswTracks");
    printf("\n");
  }

  if (Config::seedInput == cmsswSeeds && ! HasSeeds()) {
    fprintf(stderr, "Reading of CmsswSeeds requested but data not available on file.\n");
    exit(1);
  }

  if (Config::readCmsswTracks && ! HasCmsswTracks()) {
    fprintf(stderr, "Reading of CmsswTracks requested but data not available on file.\n");
    exit(1);
  }

  return f_header.f_n_events;
}

void DataFile::OpenWrite(const std::string& fname, int nev, int extra_sections)
{
  f_fp = fopen(fname.c_str(), "w");

  f_header.f_n_events = nev;

  f_header.f_extra_sections = extra_sections;

  fwrite(&f_header, sizeof(DataFileHeader), 1, f_fp);
}

int DataFile::AdvancePosToNextEvent(FILE *fp)
{
  int evsize;

  std::lock_guard<std::mutex> readlock(f_next_ev_mutex);

  fseek(fp, f_pos, SEEK_SET);
  fread(&evsize, sizeof(int), 1, fp);
  if(Config::loopOverFile) {
    // File ended, rewind back to beginning
    if(feof(fp) != 0) {
      f_pos = sizeof(DataFileHeader);
      fseek(fp, f_pos, SEEK_SET);
      fread(&evsize, sizeof(int), 1, fp);
    }
  }

  f_pos += evsize;

  return evsize;
}

void DataFile::SkipNEvents(int n_to_skip)
{
  int evsize;

  std::lock_guard<std::mutex> readlock(f_next_ev_mutex);

  while (n_to_skip-- > 0)
  {
    fseek(f_fp, f_pos, SEEK_SET);
    fread(&evsize, sizeof(int), 1, f_fp);
    f_pos += evsize;
  }
}

void DataFile::Close()
{
  fclose(f_fp);
  f_fp = 0;
  f_header = DataFileHeader();
}

void DataFile::CloseWrite(int n_written){
  if (f_header.f_n_events != n_written){
    fseek(f_fp, 0, SEEK_SET);
    f_header.f_n_events = n_written;
    fwrite(&f_header, sizeof(DataFileHeader), 1, f_fp);
  }
  Close();
}

} // end namespace mkfit
