#include "buildtestMPlex.h"

#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"

#include "MkBuilder.h"

#ifdef USE_CUDA
#include "FitterCU.h"
#include "BuilderCU.h"
#include "check_gpu_hit_structures.h"
#endif

#if defined(USE_VTUNE_PAUSE)
#include "ittnotify.h"
#endif

#include <memory>

namespace mkfit {

inline bool sortByHitsChi2(const std::pair<Track, TrackState>& cand1,
                           const std::pair<Track, TrackState>& cand2)
{
  if (cand1.first.nFoundHits() == cand2.first.nFoundHits())
    return cand1.first.chi2() < cand2.first.chi2();

  return cand1.first.nFoundHits() > cand2.first.nFoundHits();
}

inline bool sortByScore(const std::pair<Track, TrackState>& cand1,
			const std::pair<Track, TrackState>& cand2)
{
  return sortByScoreCandPair(cand1,cand2);
}

inline bool sortByPhi(const Hit& hit1, const Hit& hit2)
{
  return std::atan2(hit1.y(),hit1.x()) < std::atan2(hit2.y(),hit2.x());
}

inline bool sortByEta(const Hit& hit1, const Hit& hit2)
{
  return hit1.eta()<hit2.eta();
}

inline bool sortTracksByEta(const Track& track1, const Track& track2)
{
  return track1.momEta() < track2.momEta();
}

inline bool sortTracksByPhi(const Track& track1, const Track& track2)
{
  return track1.momPhi() < track2.momPhi();
}

struct sortTracksByPhiStruct
{
  const std::vector<std::vector<Track>>& m_track_candidates;

  sortTracksByPhiStruct(std::vector<std::vector<Track>>* track_candidates)
    : m_track_candidates( * track_candidates)
  {}

  bool operator() (const std::pair<int,int>& track1, const std::pair<int,int>& track2)
  {
    return m_track_candidates[track1.first][track1.second].posPhi() <
           m_track_candidates[track2.first][track2.second].posPhi();
  }
};

// within a layer with a "reasonable" geometry, ordering by Z is the same as eta
inline bool sortByZ(const Hit& hit1, const Hit& hit2)
{
  return hit1.z() < hit2.z();
}

//==============================================================================
// NaN and Silly track parameter check
//==============================================================================

namespace
{

  int check_nan_n_silly(TrackVec &tracks, const char* prefix)
  {
    int count = 0;
    for (auto & t : tracks)
    {
      if (t.hasSillyValues(Config::nan_n_silly_print_bad_cands_bkfit,
                           false, prefix))
      {
        ++count;
      }
    }
    return count;
  }

  void check_nan_n_silly_candiates(Event &ev)
  {
    if (Config::nan_n_silly_check_cands_every_layer)
    {
      int sc = (int) ev.nan_n_silly_per_layer_count_;
      if (sc > 0)
        printf("Nan'n'Silly: Number of silly candidates over all layers = %d\n", sc);
    }
    if (Config::nan_n_silly_check_cands_pre_bkfit)
    {
      int sc = check_nan_n_silly(ev.candidateTracks_, "Pre-bkfit silly check");
      if (sc > 0)
        printf("Nan'n'Silly: Number of silly pre-bkfit candidates = %d\n", sc);
    }
  }

  void check_nan_n_silly_bkfit(Event &ev)
  {
    if (Config::nan_n_silly_check_cands_post_bkfit)
    {
      int sc = check_nan_n_silly(ev.fitTracks_, "Post-bkfit silly check");
      if (sc > 0)
        printf("Nan'n'Silly: Number of silly post-bkfit candidates = %d\n", sc);
    }
  }

}

//==============================================================================
// runBuildTestPlexDumbCMSSW
//==============================================================================

void runBuildingTestPlexDumbCMSSW(Event& ev, MkBuilder& builder)
{
  builder.begin_event(&ev, __func__);

  if (Config::sim_val_for_cmssw) {
    builder.root_val_dumb_cmssw();
  }

  builder.end_event();
}

//==============================================================================
// runBuildTestPlexBestHit
//==============================================================================

double runBuildingTestPlexBestHit(Event& ev, MkBuilder& builder)
{
  builder.begin_event(&ev, __func__);

  builder.PrepareSeeds();

  // EventOfCandidates event_of_cands;
  builder.find_tracks_load_seeds_BH();

#ifdef USE_VTUNE_PAUSE
  __SSC_MARK(0x111);  // use this to resume Intel SDE at the same point
  __itt_resume();
#endif

#if USE_CUDA
  //check_event_of_hits_gpu(builder.get_event_of_hits());
  //check_event_of_cands_gpu(event_of_cands);
  BuilderCU builder_cu;
  builder_cu.setUpBH(builder.get_event_of_hits(), builder.get_event(),
                     event_of_cands);
#endif

  double time = dtime();

#if USE_CUDA
  builder_cu.FindTracksBestHit(event_of_cands);
  builder_cu.tearDownBH();
#else
  builder.FindTracksBestHit();
#endif

  time = dtime() - time;

#ifdef USE_VTUNE_PAUSE
  __itt_pause();
  __SSC_MARK(0x222);  // use this to pause Intel SDE at the same point
#endif

  //For best hit, the candidateTracks_ vector is the direct input to the backward fit so only need to do find_duplicates once
  if (Config::quality_val || Config::sim_val || Config::cmssw_val || Config::cmssw_export)
  {
    //Mark tracks as duplicates; if within CMSSW, remove duplicate tracks before backward fit   
    if(Config::removeDuplicates)
    {
      builder.find_duplicates(ev.candidateTracks_);
      if(Config::cmssw_export) builder.remove_duplicates(ev.candidateTracks_);
    }
  }

  // now do backwards fit... do we want to time this section?
  if (Config::backwardFit)
  {
    builder.BackwardFitBH();
  }

  if        (Config::quality_val) {
    builder.quality_val();
  } else if (Config::sim_val || Config::cmssw_val) {
    builder.root_val();
  }
  builder.end_event();
  
  // ev.print_tracks(ev.candidateTracks_, true);

  return time;
}

#if USE_CUDA
double runBuildingTestPlexBestHitGPU(Event& ev, MkBuilder& builder,
                                     BuilderCU& builder_cu)
{
  builder.begin_event(&ev, 0, __func__);

  if   (Config::seedInput == findSeeds) {builder.find_seeds();}
  else                                  {builder.map_seed_hits();} // all other simulated seeds need to have hit indices line up in LOH for seed fit

  builder.fit_seeds();

  EventOfCandidates event_of_cands;
  builder.find_tracks_load_seeds(event_of_cands);
  // Allocate event specific arrays
  builder_cu.setUpBH(builder.get_event_of_hits(), builder.get_event(),
                     event_of_cands);

  double time = dtime();

  builder_cu.FindTracksBestHit(event_of_cands);
  // Deallocate event specific arrays 
  builder_cu.tearDownBH();

  time = dtime() - time;
  if   (Config::quality_val) {
    if (!Config::silent) builder.quality_val_BH(event_of_cands);
  } else {
    builder.sim_val_BH(event_of_cands);
  }

  builder.end_event();
  
  return time;
}
#endif


//==============================================================================
// runBuildTestPlex Combinatorial: Standard TBB
//==============================================================================

double runBuildingTestPlexStandard(Event& ev, MkBuilder& builder)
{
  builder.begin_event(&ev, __func__);

  builder.PrepareSeeds();

  builder.find_tracks_load_seeds();

#ifdef USE_VTUNE_PAUSE
  __SSC_MARK(0x111);  // use this to resume Intel SDE at the same point
  __itt_resume();
#endif

  double time = dtime();

  builder.FindTracksStandard();

  time = dtime() - time;

#ifdef USE_VTUNE_PAUSE
  __itt_pause();
  __SSC_MARK(0x222);  // use this to pause Intel SDE at the same point
#endif

  check_nan_n_silly_candiates(ev);

  // first store candidate tracks
  if (Config::quality_val || Config::sim_val || Config::cmssw_val || Config::cmssw_export)
  {
    builder.quality_store_tracks(ev.candidateTracks_);
  }

  // now do backwards fit... do we want to time this section?
  if (Config::backwardFit)
  {
    builder.BackwardFit();

    check_nan_n_silly_bkfit(ev);

    if (Config::sim_val || Config::cmssw_val || Config::cmssw_export)
    {
      builder.quality_store_tracks(ev.fitTracks_);
    }
  }

  builder.handle_duplicates();

  // validation section
  if        (Config::quality_val) {
    builder.quality_val();
  } else if (Config::sim_val || Config::cmssw_val) { 
    builder.root_val();
  } else if (Config::cmssw_export) {
    builder.cmssw_export();
  }

  builder.end_event();

  // ev.print_tracks(ev.candidateTracks_, true);

  return time;
}

//==============================================================================
// runBuildTestPlex Combinatorial: CloneEngine TBB
//==============================================================================

double runBuildingTestPlexCloneEngine(Event& ev, MkBuilder& builder)
{
  builder.begin_event(&ev, __func__);

  builder.PrepareSeeds();

  builder.find_tracks_load_seeds();

#ifdef USE_VTUNE_PAUSE
  __SSC_MARK(0x111);  // use this to resume Intel SDE at the same point
  __itt_resume();
#endif
  double time = dtime();

  builder.FindTracksCloneEngine();

  time = dtime() - time;

#ifdef USE_VTUNE_PAUSE
  __itt_pause();
  __SSC_MARK(0x222);  // use this to pause Intel SDE at the same point
#endif

  check_nan_n_silly_candiates(ev);

  // first store candidate tracks
  if (Config::quality_val || Config::sim_val || Config::cmssw_val || Config::cmssw_export)
  {
    builder.quality_store_tracks(ev.candidateTracks_);
  }

  // now do backwards fit... do we want to time this section?
  if (Config::backwardFit)
  {
    builder.BackwardFit();

    check_nan_n_silly_bkfit(ev);

    if (Config::sim_val || Config::cmssw_val || Config::cmssw_export)
    {
      builder.quality_store_tracks(ev.fitTracks_);
    }
  }

  builder.handle_duplicates();

  // validation section
  if        (Config::quality_val) {
    builder.quality_val();
  } else if (Config::sim_val || Config::cmssw_val) { 
    builder.root_val();
  } else if (Config::cmssw_export) {
    builder.cmssw_export();
  }

  builder.end_event();

  // ev.print_tracks(ev.candidateTracks_, true);

  return time;
}

//==============================================================================
// runBuildTestPlex Combinatorial: Full Vector TBB
//==============================================================================

double runBuildingTestPlexFV(Event& ev, MkBuilder& builder)
{
  builder.begin_event(&ev, __func__);

  builder.PrepareSeeds();

  builder.find_tracks_load_seeds();

#ifdef USE_VTUNE_PAUSE
  __SSC_MARK(0x111);  // use this to resume Intel SDE at the same point
  __itt_resume();
#endif
  double time = dtime();

  builder.FindTracksFV();

  time = dtime() - time;

#ifdef USE_VTUNE_PAUSE
  __itt_pause();
  __SSC_MARK(0x222);  // use this to pause Intel SDE at the same point
#endif

  // first store candidate tracks
  if (Config::quality_val || Config::sim_val || Config::cmssw_val || Config::cmssw_export)
  {
    builder.quality_store_tracks(ev.candidateTracks_);
  }

  // now do backwards fit... do we want to time this section?
  if (Config::backwardFit)
  {
    builder.BackwardFit();

    if (Config::sim_val || Config::cmssw_val || Config::cmssw_export)
    {
      builder.quality_store_tracks(ev.fitTracks_);
    }
  }

  builder.handle_duplicates();

  // validation section
  if        (Config::quality_val) {
    builder.quality_val();
  } else if (Config::sim_val || Config::cmssw_val) { 
    builder.root_val();
  } else if (Config::cmssw_export) {
    builder.cmssw_export();
  }

  builder.end_event();

  // ev.print_tracks(ev.candidateTracks_, true);

  return time;
}

#if USE_CUDA
double runBuildingTestPlexCloneEngineGPU(Event& ev, EventTmp& ev_tmp, 
                                         MkBuilder& builder,
                                         BuilderCU& builder_cu,
                                         bool seed_based)
{
  EventOfCombCandidates &event_of_comb_cands = ev_tmp.m_event_of_comb_cands;
  event_of_comb_cands.Reset();

  builder.begin_event(&ev, &ev_tmp, __func__);

  if   (Config::seedInput == findSeeds) {builder.find_seeds();}
  else                                  {builder.map_seed_hits();} // all other simulated seeds need to have hit indices line up in LOH for seed fit

  builder.fit_seeds();

  builder.find_tracks_load_seeds();

#ifdef USE_VTUNE_PAUSE
  __itt_resume();
#endif

  //builder_cu.setUpFitterCE(-1 [> does not matter for now <]);
  builder_cu.allocateCE(builder.get_event_of_hits(), builder.get_event(),
                      event_of_comb_cands);
  builder_cu.setUpCE(builder.get_event_of_hits(), builder.get_event(),
                     event_of_comb_cands);

  double time = dtime();

  //builder.FindTracksCloneEngine();
  builder_cu.FindTracksCloneEngine(event_of_comb_cands, seed_based);

  time = dtime() - time;

  builder_cu.tearDownCE();
#ifdef USE_VTUNE_PAUSE
  __itt_pause();
#endif

  if (Config::quality_val) {
    if (!Config::silent) builder.quality_val_COMB();
  } else {builder.sim_val_COMB();}

  builder.end_event();

  return time;
}
#endif

//==============================================================================
// runAllBuildTestPlexBestHitGPU
//==============================================================================

#if USE_CUDA
double runAllBuildingTestPlexBestHitGPU(std::vector<Event> &events)
{

  int num_builders = events.size();
  std::vector<std::unique_ptr<MkBuilder>> builder_ptrs(num_builders);
  std::vector<EventOfCandidates> event_of_cands_vec(num_builders);
  std::vector<BuilderCU> builder_cu_vec(num_builders);

  for (int i = 0; i < builder_ptrs.size(); ++i) {
    Event &ev = events[i];
    builder_ptrs[i] = std::unique_ptr<MkBuilder> (MkBuilder::make_builder());

    MkBuilder &builder = * builder_ptrs[i].get();

    builder.begin_event(&ev, 0, __func__);

    if   (Config::seedInput == findSeeds) {builder.find_seeds();}
    else                                  {builder.map_seed_hits();} // all other simulated seeds need to have hit indices line up in LOH for seed fit

    builder.fit_seeds();

    EventOfCandidates &event_of_cands = event_of_cands_vec[i];
    builder.find_tracks_load_seeds(event_of_cands);

    BuilderCU &builder_cu = builder_cu_vec[i];
    builder_cu.setUpBH(builder.get_event_of_hits(), builder.get_event(),
                       event_of_cands);
  }

  //omp_set_num_threads(Config::numThreadsEvents);
  //std::cerr << "num threads "<< omp_get_num_threads() << std::endl;
//#pragma omp parallel for reduction(+:total_time)
  //for (int i = 0; i < builder_ptrs.size(); ++i) {
  double time = dtime();
  tbb::parallel_for(size_t(0), builder_ptrs.size(), [&](size_t i) {
    EventOfCandidates &event_of_cands = event_of_cands_vec[i];
    BuilderCU &builder_cu = builder_cu_vec[i];
    MkBuilder &builder = * builder_ptrs[i].get();

    builder_cu.FindTracksBestHit(event_of_cands);
  });
  time = dtime() - time;

  for (int i = 0; i < builder_ptrs.size(); ++i) {
    EventOfCandidates &event_of_cands = event_of_cands_vec[i];
    BuilderCU &builder_cu = builder_cu_vec[i];
    builder_cu.tearDownBH();
    MkBuilder &builder = * builder_ptrs[i].get();
    if   (!Config::sim_val && !Config::cmssw_val) {
      if (!Config::silent) builder.quality_val();
    } else if (Config::sim_val) {
      builder.root_val();
    }

    builder.end_event();
  }
  
  return time;
}
#endif
} // end namespace mkfit
