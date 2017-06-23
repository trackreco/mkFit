#include <memory>

#include "MkBuilder.h"
#include "MkBuilderEndcap.h"
#include "seedtestMPlex.h"

#include "Event.h"
#include "EventTmp.h"

#include "MkFitter.h"

//#define DEBUG
#include "Debug.h"

#include <tbb/tbb.h>

ExecutionContext g_exe_ctx;

namespace
{
  auto retcand = [](CandCloner* cloner) { g_exe_ctx.m_cloners.ReturnToPool(cloner); };
  auto retfitr = [](MkFitter*   mkfp  ) { g_exe_ctx.m_fitters.ReturnToPool(mkfp);   };
}

MkBuilder* MkBuilder::make_builder()
{
  if (Config::endcapTest) return new MkBuilderEndcap;
  else                    return new MkBuilder;
}

#ifdef DEBUG
namespace {
  void pre_prop_print(int ilay, MkFitter* mkfp) {
    std::cout << "propagate to lay=" << ilay+1 
              << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)
              << " r=" << getHypot(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1))
              << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5)
#ifdef CCSCOORD
              << " pT=" << 1./mkfp->getPar(0, 0, 3) << std::endl;
#else
              << " pT=" << getHypot(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4)) << std::endl;
#endif
  }

  void post_prop_print(int ilay, MkFitter* mkfp) {
    std::cout << "propagate to lay=" << ilay+1 
              << " arrive at x=" << mkfp->getPar(0, 1, 0) << " y=" << mkfp->getPar(0, 1, 1) << " z=" << mkfp->getPar(0, 1, 2)
              << " r=" << getHypot(mkfp->getPar(0, 1, 0), mkfp->getPar(0, 1, 1)) << std::endl;
  }

  void print_seed(const Track& seed) {
    std::cout << "MX - found seed with nHits=" << seed.nFoundHits() << " chi2=" << seed.chi2()
              << " posEta=" << seed.posEta() << " posPhi=" << seed.posPhi() << " posR=" << seed.posR()
              << " pT=" << seed.pT() << std::endl;
  }

  void print_seed2(const Track& seed) {
    std::cout << "MX - found seed with nFoundHits=" << seed.nFoundHits() << " chi2=" << seed.chi2() 
              << " x=" << seed.x() << " y=" << seed.y() << " z=" << seed.z()
              << " px=" << seed.px() << " py=" << seed.py() << " pz=" << seed.pz()
              << " pT=" << seed.pT() << std::endl;
  }

  void print_seeds(const TrackVec& seeds) {
    std::cout << "found total seeds=" << seeds.size() << std::endl;
    for (auto&& seed : seeds) {
      print_seed(seed);
    }
  }

  void print_seeds(const EventOfCandidates& event_of_cands) {
    for (int ebin = 0; ebin < Config::nEtaBin; ++ebin) {
      const EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin]; 
      for (int iseed = 0; iseed < etabin_of_candidates.m_fill_index; iseed++) {
        print_seed2(etabin_of_candidates.m_candidates[iseed]);
      }
    }
  }

  void print_seeds(const EventOfCombCandidates& event_of_comb_cands) {
    for (int ebin = 0; ebin < Config::nEtaBin; ++ebin) {
      const EtaBinOfCombCandidates &etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin]; 
      for (int iseed = 0; iseed < etabin_of_comb_candidates.m_fill_index; iseed++) {
        print_seed2(etabin_of_comb_candidates.m_candidates[iseed].front());
      }
    }
  }
}
#endif

namespace
{
  bool sortCandByHitsChi2(const Track& cand1, const Track& cand2)
  {
    if (cand1.nFoundHits() == cand2.nFoundHits())
      return cand1.chi2() < cand2.chi2();

    return cand1.nFoundHits() > cand2.nFoundHits();
  }
}

//------------------------------------------------------------------------------
// Constructor and destructor
//------------------------------------------------------------------------------

MkBuilder::MkBuilder() :
  m_event(0),
  m_event_tmp(0),
  m_event_of_hits(Config::nLayers)
{
}

MkBuilder::~MkBuilder()
{
}

//------------------------------------------------------------------------------
// Common functions
//------------------------------------------------------------------------------

void MkBuilder::begin_event(Event* ev, EventTmp* ev_tmp, const char* build_type)
{
  m_event     = ev;
  m_event_tmp = ev_tmp;

  std::vector<Track>& simtracks = m_event->simTracks_;
  // DDDD MT: debug seed fit divergence between host / mic.
  // Use this once you know seed index + set debug in MkFitter.cc, PropagationXX.cc, KalmanUtils.cc
  // Track xx = simtracks[2069];
  // simtracks.clear();
  // simtracks.push_back(xx);

  if (!Config::silent) {
    std::cout << "Building tracks with '" << build_type << "', total simtracks=" << simtracks.size() << std::endl;
  }
#ifdef DEBUG
  //unit test for eta partitioning
  for (int i = 0; i < 60; ++i)
  {
    const float eta = -1.5f + 0.05f*i;
    int b1, b2;
    int cnt = getBothEtaBins(eta, b1, b2);
    dprint("eta=" << eta << " bin=" << getEtaBin(eta) << " hb1=" << b1 << " hb2=" << b2);
  }
  //dump sim tracks
  for (int itrack = 0; itrack < simtracks.size(); ++itrack)
  {
    Track track = simtracks[itrack];
    //if (track.label() != itrack)
    //{
    //dprintf("Bad label for simtrack %d -- %d\n", itrack, track.label());
    //}
    dprint("MX - simtrack with nHits=" << track.nFoundHits() << " chi2=" << track.chi2()
              << " pT=" << track.pT() <<" phi="<< track.momPhi() <<" eta=" << track.momEta());
  }
#endif

  m_event_of_hits.Reset();

  //fill vector of hits in each layer
  tbb::parallel_for(tbb::blocked_range<int>(0, m_event->layerHits_.size()),
    [&](const tbb::blocked_range<int>& layers)
  {
    for (int ilay = layers.begin(); ilay < layers.end(); ++ilay)
    {
      m_event_of_hits.SuckInHits(m_event->layerHits_[ilay], ilay);
    }
  });

#ifdef DEBUG
  for (int itrack = 0; itrack < simtracks.size(); ++itrack)
  {
    for (int ihit = 0; ihit < simtracks[itrack].nFoundHits(); ++ihit)
    {
      dprint("track #" << itrack << " hit #" << ihit+1
	            << " hit pos=" << simtracks[itrack].hitsVector(m_event->layerHits_)[ihit].position()
              << " phi=" << simtracks[itrack].hitsVector(m_event->layerHits_)[ihit].phi()
              << " phiPart=" << getPhiPartition(simtracks[itrack].hitsVector(m_event->layerHits_)[ihit].phi()));
    }
  }
#endif

  // for (int l=0; l<m_event_of_hits.m_layers_of_hits.size(); ++l) {
  //   for (int eb=0; eb<m_event_of_hits.m_layers_of_hits[l].m_bunches_of_hits.size(); ++eb) {
  //     std::cout << "l=" << l << " eb=" << eb << " m_fill_index=" << m_event_of_hits.m_layers_of_hits[l].m_bunches_of_hits[eb].m_fill_index << " m_fill_index_old=" << m_event_of_hits.m_layers_of_hits[l].m_bunches_of_hits[eb].m_fill_index_old << std::endl;      
  //     for (int pb=0; pb<m_event_of_hits.m_layers_of_hits[l].m_bunches_of_hits[eb].m_phi_bin_infos.size(); ++pb) {
  //     	std::cout << "l=" << l << " eb=" << eb << " pb=" << pb << " first=" << m_event_of_hits.m_layers_of_hits[l].m_bunches_of_hits[eb].m_phi_bin_infos[pb].first << " second=" <<  m_event_of_hits.m_layers_of_hits[l].m_bunches_of_hits[eb].m_phi_bin_infos[pb].second << std::endl;
  //     }
  //   }
  // }

  if (! (Config::readCmsswSeeds || Config::findSeeds)) m_event->seedTracks_ = m_event->simTracks_; // make seed tracks == simtracks if not using "realistic" seeding
}

void MkBuilder::end_event()
{
  m_event = 0;
}

//------------------------------------------------------------------------------
// Seeding functions: finding and fitting
//------------------------------------------------------------------------------

int MkBuilder::find_seeds()
{
#ifdef DEBUG
  bool debug(false);
#endif
  TripletIdxConVec seed_idcs;

  double time = dtime();
  findSeedsByRoadSearch(seed_idcs,m_event_of_hits.m_layers_of_hits,m_event->layerHits_[1].size(),m_event);
  time = dtime() - time;

  // use this to initialize tracks
  const Hit * lay0hits = m_event_of_hits.m_layers_of_hits[0].m_hits;
  const Hit * lay1hits = m_event_of_hits.m_layers_of_hits[1].m_hits;
  const Hit * lay2hits = m_event_of_hits.m_layers_of_hits[2].m_hits;

  // make seed tracks
  TrackVec & seedtracks = m_event->seedTracks_;
  seedtracks.resize(seed_idcs.size());
  for (int iseed = 0; iseed < seedtracks.size(); iseed++)
  {
    auto & seedtrack = seedtracks[iseed];
    seedtrack.setLabel(iseed);

    // use to set charge
    const Hit & hit0 = lay0hits[seed_idcs[iseed][0]];
    const Hit & hit1 = lay1hits[seed_idcs[iseed][1]];
    const Hit & hit2 = lay2hits[seed_idcs[iseed][2]];

    seedtrack.setCharge(calculateCharge(hit0,hit1,hit2));

    for (int ihit = 0; ihit < Config::nlayers_per_seed; ihit++)
    {
      seedtrack.addHitIdx(seed_idcs[iseed][ihit],0.0f);
    }

    for (int ihit = Config::nlayers_per_seed; ihit < Config::nLayers; ihit++)
    {
      seedtrack.setHitIdx(ihit,-1);
    }
    
    dprint("iseed: " << iseed << " mcids: " << hit0.mcTrackID(m_event->simHitsInfo_) << " " <<
	   hit1.mcTrackID(m_event->simHitsInfo_) << " " << hit1.mcTrackID(m_event->simHitsInfo_));
  }
  return time;
}

void MkBuilder::fit_seeds()
{
  g_exe_ctx.populate(Config::numThreadsFinder);
  TrackVec& seedtracks = m_event->seedTracks_;

  int theEnd = seedtracks.size();
  int count = (theEnd + NN - 1)/NN;

  tbb::parallel_for(tbb::blocked_range<int>(0, count, std::max(1, Config::numSeedsPerTask/NN)),
    [&](const tbb::blocked_range<int>& i) {

      std::unique_ptr<MkFitter, decltype(retfitr)> mkfp(g_exe_ctx.m_fitters.GetFromPool(), retfitr);
      for (int it = i.begin(); it < i.end(); ++it)
      {
        fit_one_seed_set(seedtracks, it*NN, std::min((it+1)*NN, theEnd), mkfp.get());
      }
    }
  );

  //ok now, we should have all seeds fitted in recseeds
  dcall(print_seeds(seedtracks));
}

inline void MkBuilder::fit_one_seed_set(TrackVec& seedtracks, int itrack, int end, MkFitter *mkfp)
{
  mkfp->SetNhits(Config::nlayers_per_seed); //just to be sure (is this needed?)
  mkfp->InputTracksAndHits(seedtracks, m_event_of_hits.m_layers_of_hits, itrack, end);
  if (Config::cf_seeding) mkfp->ConformalFitTracks(false, itrack, end);
  if (Config::readCmsswSeeds==false) mkfp->FitTracks(end - itrack, m_event);

  const int ilay = Config::nlayers_per_seed; // layer 4

  dcall(pre_prop_print(ilay, mkfp));
  mkfp->PropagateTracksToR(m_event->geom_.Radius(ilay), end - itrack);
  dcall(post_prop_print(ilay, mkfp));

  mkfp->OutputFittedTracksAndHitIdx(m_event->seedTracks_, itrack, end, true);
}

//------------------------------------------------------------------------------
// Common functions for validation
//------------------------------------------------------------------------------

////////////////////////////////
// Outline of map/remap logic //
////////////////////////////////
/* 
All built candidate tracks have all hit indices pointing to m_event_of_hits.m_layers_of_hits[layer].m_hits (LOH)
MC seeds (both CMSSW and toyMC) have seed hit indices pointing to global HitVec m_event->layerHits_[layer] (GLH)
"Real" seeds have all seed hit indices pointing to LOH.
So.. to have universal seed fitting function --> have seed hits point to LOH no matter their origin.
This means that all MC seeds must be "mapped" from GLH to LOH: map_seed_hits().
Now InputTracksAndHits() for seed fit will use LOH instead of GLH.
The output tracks of the seed fitting are now stored in m_event->seedTracks_.

Then building proceeds as normal, using m_event->seedTracks_ as input no matter the choice of seeds. 

For the validation, we can reuse the TrackExtra setMCTrackIDInfo() with a few tricks.
Since setMCTrackIDInfo by necessity uses GLH, we then need ALL track collections (seed, candidate, fit) to their hits point back to GLH.
There are also two validation options: w/ or w/o ROOT.

W/ ROOT uses the TTreValidation class which needs seedTracks_, candidateTracks_, and fitTracks_ all stored in m_event.
The fitTracks_ collection for now is just a copy of candidateTracks_ (eventually may have cuts and things that affect which tracks to fit).
So... need to "remap" seedTracks_ hits from LOH to GLH with remap_seed_hits().
And also copy in tracks from EtaBin* to candidateTracks_, and then remap hits from LOH to GLH with quality_store_tracks() and remap_cand_hits().
W/ ROOT uses root_val_BH for BH, and root_val_COMB() for non-BH.

W/O ROOT is a bit simpler... as we only need to do the copy out tracks from EtaBin* and then remap just candidateTracks_.
This uses quality_output_COMB() or quality_output_BH()

N.B.1 Since fittestMPlex at the moment is not "end-to-end" with candidate tracks, we can still use the GLH version of InputTracksAndHits()

N.B.2 Since we inflate LOH by 2% more than GLH, hit indices in building only go to GLH, so all loops are sized to GLH.
*/

void MkBuilder::map_seed_hits()
{ // map seed hit indices from global m_event->layerHits_[i] to hit indices in structure m_event_of_hits.m_layers_of_hits[i].m_hits
  HitIDVec seedLayersHitMap(m_event->simHitsInfo_.size());
  for (int ilayer = 0; ilayer < Config::nlayers_per_seed; ++ilayer) {
    const auto & lof_m_hits = m_event_of_hits.m_layers_of_hits[ilayer].m_hits;
    const auto   size = m_event->layerHits_[ilayer].size();
    for (int index = 0; index < size; ++index) { 
      seedLayersHitMap[lof_m_hits[index].mcHitID()] = HitID(ilayer, index);
    }
  }
  for (int ilayer = 0; ilayer < Config::nlayers_per_seed; ++ilayer) {
    const auto & global_hit_vec = m_event->layerHits_[ilayer];
    const auto   size = m_event->layerHits_[ilayer].size();
    for (auto&& track : m_event->seedTracks_) {
      auto hitidx = track.getHitIdx(ilayer);
      if ((hitidx>=0) && (hitidx<size)) 
      {
	track.setHitIdx(ilayer, seedLayersHitMap[global_hit_vec[hitidx].mcHitID()].index);
      }
    }
  }
}

void MkBuilder::remap_seed_hits()
{ // map seed hit indices from hit indices in structure m_event_of_hits.m_layers_of_hits[i].m_hits to global m_event->layerHits_[i]
  HitIDVec seedLayersHitMap(m_event->simHitsInfo_.size());
  for (int ilayer = 0; ilayer < Config::nlayers_per_seed; ++ilayer) {
    const auto & global_hit_vec = m_event->layerHits_[ilayer];
    const auto   size = global_hit_vec.size();
    for (int index = 0; index < size; ++index) {
      seedLayersHitMap[global_hit_vec[index].mcHitID()] = HitID(ilayer, index);
    }
  }
  for (int ilayer = 0; ilayer < Config::nlayers_per_seed; ++ilayer) {
    const auto & lof_m_hits = m_event_of_hits.m_layers_of_hits[ilayer].m_hits;
    const auto   size = m_event->layerHits_[ilayer].size();
    for (auto&& track : m_event->seedTracks_) {
      auto hitidx = track.getHitIdx(ilayer);
      if ((hitidx>=0) && (hitidx<size))
      { 
	track.setHitIdx(ilayer, seedLayersHitMap[lof_m_hits[hitidx].mcHitID()].index);
      }
    }
  }
}

void MkBuilder::remap_cand_hits()
{ // map cand hit indices from hit indices in structure m_event_of_hits.m_layers_of_hits[i].m_hits to global m_event->layerHits_[i]
  HitIDVec candLayersHitMap(m_event->simHitsInfo_.size());
  for (int ilayer = 0; ilayer < Config::nLayers; ++ilayer) {
    const auto & global_hit_vec = m_event->layerHits_[ilayer];
    for (int index = 0; index < global_hit_vec.size(); ++index) {
      candLayersHitMap[global_hit_vec[index].mcHitID()] = HitID(ilayer, index);
    }
  }
  for (int ilayer = 0; ilayer < Config::nLayers; ++ilayer) {
    const auto & lof_m_hits = m_event_of_hits.m_layers_of_hits[ilayer].m_hits;
    const auto   size = m_event->layerHits_[ilayer].size();
    for (auto&& track : m_event->candidateTracks_) {
      auto hitidx = track.getHitIdx(ilayer);
      if ((hitidx>=0) && (hitidx<size)) 
      {	
	track.setHitIdx(ilayer, candLayersHitMap[lof_m_hits[hitidx].mcHitID()].index);
      }
    }
  }
}

void MkBuilder::align_simtracks()
{
  if (Config::readCmsswSeeds && Config::endcapTest) 
  {
    for (int itrack = 0; itrack < m_event->simTracks_.size(); itrack++)
    {
      m_event->simTracks_[itrack].setLabel(itrack);
    }
  }
}

//------------------------------------------------------------------------------
// Non-ROOT validation
//------------------------------------------------------------------------------

void MkBuilder::quality_output_BH(const EventOfCandidates& event_of_cands)
{
  quality_reset();

  quality_store_tracks_BH(event_of_cands);
  
  remap_cand_hits();

  align_simtracks();

  for (int itrack = 0; itrack < m_event->candidateTracks_.size(); itrack++)
  {
    quality_process(m_event->candidateTracks_[itrack]);
  }
  
  quality_print();
}

void MkBuilder::quality_output_COMB()
{
  quality_reset();

  quality_store_tracks_COMB();

  remap_cand_hits();

  align_simtracks();

  for (int iseed = 0; iseed < m_event->candidateTracks_.size(); iseed++)
  {
    quality_process(m_event->candidateTracks_[iseed]);
  }

  quality_print();
}

void MkBuilder::quality_reset()
{
  m_cnt = m_cnt1 = m_cnt2 = m_cnt_8 = m_cnt1_8 = m_cnt2_8 = m_cnt_nomc = 0;
}

void MkBuilder::quality_store_tracks_BH(const EventOfCandidates& event_of_cands)
{
  for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
  {
    const EtaBinOfCandidates &etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin]; 
      
    for (int itrack = 0; itrack < etabin_of_candidates.m_fill_index; itrack++)
    {
      m_event->candidateTracks_.push_back(etabin_of_candidates.m_candidates[itrack]);
    }
  }
}

void MkBuilder::quality_store_tracks_COMB()
{
  for (int ebin = 0; ebin < Config::nEtaBin; ++ebin)
  {
    const EtaBinOfCombCandidates &etabin_of_comb_candidates = m_event_tmp->m_event_of_comb_cands.m_etabins_of_comb_candidates[ebin]; 
    
    for (int iseed = 0; iseed < etabin_of_comb_candidates.m_fill_index; iseed++)
    {
      // take the first one!
      if ( ! etabin_of_comb_candidates.m_candidates[iseed].empty())
      {
     	m_event->candidateTracks_.push_back(etabin_of_comb_candidates.m_candidates[iseed].front());
      }
    }
  }
}

void MkBuilder::quality_process(Track &tkcand)
{
  TrackExtra extra(tkcand.label());
  extra.setMCTrackIDInfo(tkcand, m_event->layerHits_, m_event->simHitsInfo_);
  int mctrk = extra.mcTrackID();

  float pt    = tkcand.pT();
  float ptmc = 0., pr = 0., nfoundmc = 0., chi2mc = 0.;

  if (mctrk < 0 || mctrk >= Config::nTracks)
  {
    ++m_cnt_nomc;
    // std::cout << "XX bad track idx " << mctrk << "\n";
  } else {

    ptmc  = m_event->simTracks_[mctrk].pT() ;
    pr    = pt / ptmc;
    nfoundmc = m_event->simTracks_[mctrk].nFoundHits();
    chi2mc = m_event->simTracks_[mctrk].chi2();//this is actually the number of reco hits in cmssw

    ++m_cnt;
    if (pr > 0.9 && pr < 1.1) ++m_cnt1;
    if (pr > 0.8 && pr < 1.2) ++m_cnt2;

    if (tkcand.nFoundHits() >= 0.8f*nfoundmc)
      {
	++m_cnt_8;
	if (pr > 0.9 && pr < 1.1) ++m_cnt1_8;
	if (pr > 0.8 && pr < 1.2) ++m_cnt2_8;
      }
  }

#if defined(DEBUG) || defined(PRINTOUTS_FOR_PLOTS)
  if (!Config::silent) {
    std::lock_guard<std::mutex> printlock(Event::printmutex);
    std::cout << "MX - found track with nFoundHits=" << tkcand.nFoundHits() << " chi2=" << tkcand.chi2() << " pT=" << pt <<" pTmc="<< ptmc << " nfoundmc=" << nfoundmc << " chi2mc=" << chi2mc <<" lab="<< tkcand.label() <<std::endl;
  }
#endif
  // DDDD MT: debug seed fit divergence between host / mic.
  // Use this to compare track quality.
  // printf("DDDD N_h=%d lab=%d\n", tkcand.nFoundHits(), tkcand.label());
}

void MkBuilder::quality_print()
{
  if (!Config::silent) {
    std::lock_guard<std::mutex> printlock(Event::printmutex);
    std::cout << "found tracks=" << m_cnt   << "  in pT 10%=" << m_cnt1   << "  in pT 20%=" << m_cnt2   << "     no_mc_assoc="<< m_cnt_nomc <<std::endl;
    std::cout << "  nH >= 80% =" << m_cnt_8 << "  in pT 10%=" << m_cnt1_8 << "  in pT 20%=" << m_cnt2_8 << std::endl;
  }
}

//------------------------------------------------------------------------------
// Root validation
//------------------------------------------------------------------------------

void MkBuilder::root_val_BH(const EventOfCandidates& event_of_cands)
{
  // remap seed tracks
  remap_seed_hits();

  // get the tracks ready for validation
  quality_store_tracks_BH(event_of_cands);
  remap_cand_hits();
  m_event->fitTracks_ = m_event->candidateTracks_; // fixme: hack for now. eventually fitting will be including end-to-end
  align_simtracks();
  init_track_extras();

  m_event->Validate();
}

void MkBuilder::root_val_COMB()
{
  remap_seed_hits(); // prepare seed tracks for validation

  // get the tracks ready for validation
  quality_store_tracks_COMB();
  remap_cand_hits();
  m_event->fitTracks_ = m_event->candidateTracks_; // fixme: hack for now. eventually fitting will be including end-to-end
  align_simtracks();
  init_track_extras();

  m_event->Validate();
}

void MkBuilder::init_track_extras()
{
  TrackVec      & seedtracks      = m_event->seedTracks_; 
  TrackExtraVec & seedtrackextras = m_event->seedTracksExtra_;
  for (int i = 0; i < seedtracks.size(); i++)
  {
    seedtrackextras.emplace_back(seedtracks[i].label());
  }
  m_event->validation_.alignTrackExtra(seedtracks,seedtrackextras);

  TrackVec      & candidatetracks      = m_event->candidateTracks_;
  TrackExtraVec & candidatetrackextras = m_event->candidateTracksExtra_;
  for (int i = 0; i < candidatetracks.size(); i++)
  {
    candidatetrackextras.emplace_back(candidatetracks[i].label());
  }
  m_event->validation_.alignTrackExtra(candidatetracks,candidatetrackextras);
  
  TrackVec      & fittracks      = m_event->fitTracks_;
  TrackExtraVec & fittrackextras = m_event->fitTracksExtra_;
  for (int i = 0; i < fittracks.size(); i++)
  {
    fittrackextras.emplace_back(fittracks[i].label());
  }
  m_event->validation_.alignTrackExtra(fittracks,fittrackextras);
}

//------------------------------------------------------------------------------
// FindTracksBestHit
//------------------------------------------------------------------------------

void MkBuilder::find_tracks_load_seeds(EventOfCandidates& event_of_cands)
{
  // partition recseeds into eta bins
  for (int iseed = 0; iseed < m_event->seedTracks_.size(); ++iseed)
  {
    //if (m_event->seedTracks_[iseed].label() != iseed)
    //{
    //printf("Bad label for recseed %d -- %d\n", iseed, m_event->seedTracks_[iseed].label());
    //}

    event_of_cands.InsertCandidate(m_event->seedTracks_[iseed]);
  }

  //dump seeds
  dcall(print_seeds(event_of_cands));
}

void MkBuilder::FindTracksBestHit(EventOfCandidates& event_of_cands)
{
  g_exe_ctx.populate(Config::numThreadsFinder);

  tbb::parallel_for(tbb::blocked_range<int>(0, Config::nEtaBin),
    [&](const tbb::blocked_range<int>& ebins)
  {
    for (int ebin = ebins.begin(); ebin != ebins.end(); ++ebin) {
      EtaBinOfCandidates& etabin_of_candidates = event_of_cands.m_etabins_of_candidates[ebin];

      tbb::parallel_for(tbb::blocked_range<int>(0,etabin_of_candidates.m_fill_index,Config::numSeedsPerTask), 
        [&](const tbb::blocked_range<int>& tracks)
      {
        std::unique_ptr<MkFitter, decltype(retfitr)> mkfp(g_exe_ctx.m_fitters.GetFromPool(), retfitr);
        
        for (int itrack = tracks.begin(); itrack < tracks.end(); itrack += NN) {
          int end = std::min(itrack + NN, tracks.end());

          dprint(std::endl << "processing track=" << itrack << " etabin=" << ebin << " findex=" << etabin_of_candidates.m_fill_index);

          mkfp->SetNhits(3);//just to be sure (is this needed?)
          mkfp->InputTracksAndHitIdx(etabin_of_candidates.m_candidates, itrack, end, true);

          //ok now we start looping over layers
          //loop over layers, starting from after the seed
          //consider inverting loop order and make layer outer, need to trade off hit prefetching with copy-out of candidates
          for (int ilay = Config::nlayers_per_seed; ilay < Config::nLayers; ++ilay)
          {
            LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[ilay];

            // XXX This should actually be done in some other thread for the next layer while
            // this thread is crunching the current one.
            // For now it's done in MkFitter::AddBestHit(), two loops before the data is needed.
            // for (int i = 0; i < bunch_of_hits.m_fill_index; ++i)
            // {
            //   _mm_prefetch((char*) & bunch_of_hits.m_hits[i], _MM_HINT_T1);
            // }

            mkfp->SelectHitIndices(layer_of_hits, end - itrack);

// #ifdef PRINTOUTS_FOR_PLOTS
// 	     std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
// #endif

            //make candidates with best hit
            dprint("make new candidates");
            mkfp->AddBestHit(layer_of_hits, end - itrack);
            mkfp->SetNhits(ilay + 1);  //here again assuming one hit per layer (is this needed?)

            //propagate to layer
            if (ilay + 1 < Config::nLayers)
            {
              dcall(pre_prop_print(ilay, mkfp.get()));
              mkfp->PropagateTracksToR(m_event->geom_.Radius(ilay+1), end - itrack);
              dcall(post_prop_print(ilay, mkfp.get()));
            }

          } // end of layer loop
          mkfp->OutputFittedTracksAndHitIdx(etabin_of_candidates.m_candidates, itrack, end, true);
        }
      }); // end of seed loop
    }
  }); //end of parallel section over seeds
}

//------------------------------------------------------------------------------
// FindTracksCombinatorial: Standard TBB and CloneEngine TBB 
//------------------------------------------------------------------------------

void MkBuilder::find_tracks_load_seeds()
{
  EventOfCombCandidates &event_of_comb_cands = m_event_tmp->m_event_of_comb_cands;

  for (int iseed = 0; iseed < m_event->seedTracks_.size(); ++iseed)
  {
    //if (m_event->seedTracks_[iseed].label() != iseed)
    //{
    //printf("Bad label for recseed %d -- %d\n", iseed, m_event->seedTracks_[iseed].label());
    //}
    event_of_comb_cands.InsertSeed(m_event->seedTracks_[iseed]);
  }

  //dump seeds
  dcall(print_seeds(event_of_comb_cands));
}

//------------------------------------------------------------------------------
// FindTracksCombinatorial: Standard TBB
//------------------------------------------------------------------------------

void MkBuilder::FindTracksStandard()
{
  g_exe_ctx.populate(Config::numThreadsFinder);
  EventOfCombCandidates &event_of_comb_cands = m_event_tmp->m_event_of_comb_cands;

  tbb::parallel_for(tbb::blocked_range<int>(0, Config::nEtaBin),
    [&](const tbb::blocked_range<int>& ebins)
  {
    // loop over eta bins
    for (int ebin = ebins.begin(); ebin != ebins.end(); ++ebin) {
      EtaBinOfCombCandidates& etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin];
      
      int adaptiveSPT = Config::nEtaBin*etabin_of_comb_candidates.m_fill_index/Config::numThreadsFinder/2 + 1;
      dprint("adaptiveSPT " << adaptiveSPT << " fill " << etabin_of_comb_candidates.m_fill_index);
      // loop over seeds
      tbb::parallel_for(tbb::blocked_range<int>(0, etabin_of_comb_candidates.m_fill_index, std::min(Config::numSeedsPerTask, adaptiveSPT)), 
	[&](const tbb::blocked_range<int>& seeds)
      {
	const int start_seed = seeds.begin();
	const int end_seed   = seeds.end();
	const int nseeds     = end_seed - start_seed;

	//ok now we start looping over layers
	//loop over layers, starting from after the seed
	for (int ilay = Config::nlayers_per_seed; ilay < Config::nLayers; ++ilay)
	{
	  LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[ilay];
	
	  dprint("processing lay=" << ilay+1);
	
	  // prepare unrolled vector to loop over
	  std::vector<std::pair<int,int> > seed_cand_idx;
	
	  for (int iseed = start_seed; iseed < end_seed; ++iseed)
	  {
	    std::vector<Track> &scands = etabin_of_comb_candidates.m_candidates[iseed];
	    for (int ic = 0; ic < scands.size(); ++ic)
	    {
	      if (scands[ic].getLastHitIdx() >= -1)
	      {
		seed_cand_idx.push_back(std::pair<int,int>(iseed,ic));
	      }
	    }
	  }
	  int theEndCand = seed_cand_idx.size();

	  if (theEndCand == 0) continue;

	  std::vector<std::vector<Track>> tmp_candidates(nseeds);
	  for (int iseed = 0; iseed < tmp_candidates.size(); ++iseed)
	  {
	    // XXXX MT: Tried adding 25 to reserve below as I was seeing some
	    // time spent in push_back ... but it didn't really help.
	    // We need to optimize this by throwing away and replacing the worst
	    // candidate once a better one arrives. This will also avoid sorting.
	    tmp_candidates[iseed].reserve(2*Config::maxCandsPerSeed);//factor 2 seems reasonable to start with
	  }

	  //vectorized loop
	  for (int itrack = 0; itrack < theEndCand; itrack += NN)
	  {
	    int end = std::min(itrack + NN, theEndCand);
	  
	    dprint("processing track=" << itrack);
	  
	    std::unique_ptr<MkFitter,   decltype(retfitr)> mkfp  (g_exe_ctx.m_fitters.GetFromPool(), retfitr);	
	    
	    mkfp->SetNhits(ilay);//here again assuming one hit per layer
	  
	    //fixme find a way to deal only with the candidates needed in this thread
	    mkfp->InputTracksAndHitIdx(etabin_of_comb_candidates.m_candidates,
				       seed_cand_idx, itrack, end,
				       ilay == Config::nlayers_per_seed);

	    //propagate to layer
	    if (ilay > Config::nlayers_per_seed)
	    {
	      dcall(pre_prop_print(ilay, mkfp));
	      mkfp->PropagateTracksToR(m_event->geom_.Radius(ilay), end - itrack);
	      dcall(post_prop_print(ilay, mkfp));
	    }

	    dprint("now get hit range");
	    mkfp->SelectHitIndices(layer_of_hits, end - itrack);
	  
	  //#ifdef PRINTOUTS_FOR_PLOTS
	  //std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
	  //#endif
	    
	    dprint("make new candidates");
	    mkfp->FindCandidates(layer_of_hits, tmp_candidates, start_seed, end - itrack);
	  
	  } //end of vectorized loop

	  // clean exceeding candidates per seed
	  // FIXME: is there a reason why these are not vectorized????
	  for (int is = 0; is < tmp_candidates.size(); ++is)
	  {
	    dprint("dump seed n " << is << " with input candidates=" << tmp_candidates[is].size());
	    std::sort(tmp_candidates[is].begin(), tmp_candidates[is].end(), sortCandByHitsChi2);
	  
	    if (tmp_candidates[is].size() > Config::maxCandsPerSeed)
	    {
	      dprint("erase extra candidates" << " tmp_candidates[is].size()=" << tmp_candidates[is].size()
		     << " Config::maxCandsPerSeed=" << Config::maxCandsPerSeed);
	      tmp_candidates[is].erase(tmp_candidates[is].begin() + Config::maxCandsPerSeed,
				       tmp_candidates[is].end());
	    }
	    dprint("dump seed n " << is << " with output candidates=" << tmp_candidates[is].size());
	  }
	  //now swap with input candidates
	  for (int is = 0; is < tmp_candidates.size(); ++is)
	  {
	    if (tmp_candidates[is].size() > 0)
	    {
	      // Copy the best -2 cands back to the current list.
	      int num_hits = tmp_candidates[is].size();
	    
	      if (num_hits < Config::maxCandsPerSeed)
	      {
		std::vector<Track> &ov = etabin_of_comb_candidates.m_candidates[start_seed+is];
		const int max_m2 = ov.size();
	      
		int cur_m2 = 0;
		while (cur_m2 < max_m2 && ov[cur_m2].getLastHitIdx() != -2) ++cur_m2;
		while (cur_m2 < max_m2 && num_hits < Config::maxCandsPerSeed)
	        {
		  tmp_candidates[is].push_back( ov[cur_m2++] );
		  ++num_hits;
		}
	      }

	      etabin_of_comb_candidates.m_candidates[start_seed+is].swap(tmp_candidates[is]);
	      tmp_candidates[is].clear();
	    }
	  }
	  
	} // end of layer loop
	
	// final sorting
	for (int iseed = start_seed; iseed < end_seed; ++iseed)
	{
	  std::vector<Track>& finalcands = etabin_of_comb_candidates.m_candidates[iseed];
	  if (finalcands.size() == 0) continue;
	  std::sort(finalcands.begin(), finalcands.end(), sortCandByHitsChi2);
	}
      }); // end parallel-for loop over chunk of seeds within etabin
    } // end loop over over eta bins within block
  }); // end of loop over block of eta bins
}

//------------------------------------------------------------------------------
// FindTracksCombinatorial: CloneEngine TBB
//------------------------------------------------------------------------------

void MkBuilder::FindTracksCloneEngine()
{
  g_exe_ctx.populate(Config::numThreadsFinder);
  EventOfCombCandidates &event_of_comb_cands = m_event_tmp->m_event_of_comb_cands;

  tbb::parallel_for(tbb::blocked_range<int>(0, Config::nEtaBin),
    [&](const tbb::blocked_range<int>& ebins)
  {
    for (int ebin = ebins.begin(); ebin != ebins.end(); ++ebin) {
      EtaBinOfCombCandidates& etabin_of_comb_candidates = event_of_comb_cands.m_etabins_of_comb_candidates[ebin];

      int adaptiveSPT = Config::nEtaBin*etabin_of_comb_candidates.m_fill_index/Config::numThreadsFinder/2 + 1;
      dprint("adaptiveSPT " << adaptiveSPT << " fill " << etabin_of_comb_candidates.m_fill_index);
      tbb::parallel_for(tbb::blocked_range<int>(0, etabin_of_comb_candidates.m_fill_index, std::min(Config::numSeedsPerTask, adaptiveSPT)), 
        [&](const tbb::blocked_range<int>& seeds)
      {
        std::unique_ptr<CandCloner, decltype(retcand)> cloner(g_exe_ctx.m_cloners.GetFromPool(), retcand);
        std::unique_ptr<MkFitter,   decltype(retfitr)> mkfp  (g_exe_ctx.m_fitters.GetFromPool(), retfitr);

        // loop over layers
        find_tracks_in_layers(etabin_of_comb_candidates, *cloner, mkfp.get(), seeds.begin(), seeds.end(), ebin);
      });
    }
  });
}

void MkBuilder::find_tracks_in_layers(EtaBinOfCombCandidates &etabin_of_comb_candidates, CandCloner &cloner,
                                      MkFitter *mkfp, int start_seed, int end_seed, int ebin)
{
  auto n_seeds = end_seed - start_seed;

  std::vector<std::pair<int,int>> seed_cand_idx;
  seed_cand_idx.reserve(n_seeds * Config::maxCandsPerSeed);

  cloner.begin_eta_bin(&etabin_of_comb_candidates, start_seed, n_seeds);

  //loop over layers, starting from after the seeD
  for (int ilay = Config::nlayers_per_seed; ilay <= Config::nLayers; ++ilay)
  {
    dprint("processing lay=" << ilay+1);

    //prepare unrolled vector to loop over
    for (int iseed = start_seed; iseed != end_seed; ++iseed)
    {
      std::vector<Track> &scands = etabin_of_comb_candidates.m_candidates[iseed];
      for (int ic = 0; ic < scands.size(); ++ic)
      {
        if (scands[ic].getLastHitIdx() >= -1)
        {
          seed_cand_idx.push_back(std::pair<int,int>(iseed,ic));
        }
      }
    }
    const int theEndCand = seed_cand_idx.size();

    // don't bother messing with the clone engine if there are no candidates
    // (actually it crashes, so this protection is needed)
    // XXXX MT ??? How does this happen ???
    if (theEndCand == 0) continue;

    if (ilay < Config::nLayers)
    {
      cloner.begin_layer(ilay);
    }

    //vectorized loop
    for (int itrack = 0; itrack < theEndCand; itrack += NN)
    {
      const int end = std::min(itrack + NN, theEndCand);

#ifdef DEBUG
      dprint("processing track=" << itrack);
      dprintf("FTCE: start_seed=%d, n_seeds=%d, theEndCand=%d\n"
              "      itrack=%d, end=%d, nn=%d, end_eq_tec=%d\n",
              start_seed, n_seeds, theEndCand,
              itrack, end, end-itrack, end == theEndCand);
      dprintf("      ");
      for (int i=itrack; i < end; ++i) dprintf("%d,%d  ", seed_cand_idx[i].first, seed_cand_idx[i].second);
      dprintf("\n");
#endif

      // mkfp->SetNhits(ilay == Config::nlayers_per_seed ? ilay : ilay + 1);
      mkfp->SetNhits(ilay);
      mkfp->InputTracksAndHitIdx(etabin_of_comb_candidates.m_candidates,
                                 seed_cand_idx, itrack, end,
                                 true);

#ifdef DEBUG
      for (int i=itrack; i < end; ++i)
        dprintf("  track %d, idx %d is from seed %d\n", i, i - itrack, mkfp->Label(i - itrack,0,0));
      dprintf("\n");
#endif

      if (ilay > Config::nlayers_per_seed)
      {
        LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[ilay - 1];
        mkfp->UpdateWithLastHit(layer_of_hits, end - itrack);
        // TODO: Check that: tmp arguments while porting? 
        //mkfp->UpdateWithLastHit(layer_of_hits, end - itrack, cuFitter,
            //event_of_hits_cu.m_layers_of_hits[ilay-1]);

        if (ilay < Config::nLayers)
        {
          // Propagate to this layer
          mkfp->PropagateTracksToR(m_event->geom_.Radius(ilay), end - itrack);

          // copy_out the propagated track params, errors only (hit-idcs and chi2 already updated)
          mkfp->CopyOutParErr(etabin_of_comb_candidates.m_candidates,
              end - itrack, true);
        } 
        else {
          // copy_out the updated track params, errors only (hit-idcs and chi2 already updated)
          mkfp->CopyOutParErr(etabin_of_comb_candidates.m_candidates,
              end - itrack, false);
          continue;
        }
      }

      // if (ilay == Config::nLayers)
      // {
      //   break;
      // }

      dprint("now get hit range");

      LayerOfHits &layer_of_hits = m_event_of_hits.m_layers_of_hits[ilay];

      mkfp->SelectHitIndices(layer_of_hits, end - itrack);

      //#ifdef PRINTOUTS_FOR_PLOTS
      //std::cout << "MX number of hits in window in layer " << ilay << " is " <<  mkfp->getXHitEnd(0, 0, 0)-mkfp->getXHitBegin(0, 0, 0) << std::endl;
      //#endif

      dprint("make new candidates");
      cloner.begin_iteration();

      mkfp->FindCandidatesMinimizeCopy(layer_of_hits, cloner, start_seed, end - itrack);

      cloner.end_iteration();
    } //end of vectorized loop

    if (ilay < Config::nLayers)
    {
      cloner.end_layer();
    }
    seed_cand_idx.clear();

  } // end of layer loop

  cloner.end_eta_bin();

  // final sorting
  for (int iseed = start_seed; iseed < end_seed; ++iseed)
  {
    std::vector<Track>& finalcands = etabin_of_comb_candidates.m_candidates[iseed];
    if (finalcands.size() == 0) continue;
    std::sort(finalcands.begin(), finalcands.end(), sortCandByHitsChi2);
  }
}
