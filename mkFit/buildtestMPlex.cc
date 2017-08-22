#include "buildtestMPlex.h"

#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"
#include "BinInfoUtils.h"

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

inline bool sortByHitsChi2(const std::pair<Track, TrackState>& cand1,
                           const std::pair<Track, TrackState>& cand2)
{
  if (cand1.first.nFoundHits() == cand2.first.nFoundHits())
    return cand1.first.chi2() < cand2.first.chi2();

  return cand1.first.nFoundHits() > cand2.first.nFoundHits();
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
// runBuildTestPlexBestHit
//==============================================================================

double runBuildingTestPlexBestHit(Event& ev, MkBuilder& builder)
{
  builder.begin_event(&ev, 0, __func__);

  if   (Config::findSeeds) {builder.find_seeds();}
  else                     {builder.map_seed_hits();} // all other simulated seeds need to have hit indices line up in LOH for seed fit

  builder.fit_seeds();

  EventOfCandidates event_of_cands;
  builder.find_tracks_load_seeds(event_of_cands);

#ifdef USE_VTUNE_PAUSE
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
  builder.FindTracksBestHit(event_of_cands);
#endif

  time = dtime() - time;

#ifdef USE_VTUNE_PAUSE
  __itt_pause();
#endif
  
  if   (!Config::normal_val) {
    if (!Config::silent) builder.quality_output_BH(event_of_cands);
  } else {
    builder.root_val_BH(event_of_cands);
  }

  builder.end_event();
  
  return time;
}

#if USE_CUDA
double runBuildingTestPlexBestHitGPU(Event& ev, MkBuilder& builder,
                                     BuilderCU& builder_cu)
{
  builder.begin_event(&ev, 0, __func__);

  if   (Config::findSeeds) {builder.find_seeds();}
  else                     {builder.map_seed_hits();} // all other simulated seeds need to have hit indices line up in LOH for seed fit

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
  if   (!Config::normal_val) {
    if (!Config::silent) builder.quality_output_BH(event_of_cands);
  } else {
    builder.root_val_BH(event_of_cands);
  }

  builder.end_event();
  
  return time;
}
#endif


//==============================================================================
// runBuildTestPlex Combinatorial: Standard TBB
//==============================================================================

double runBuildingTestPlexStandard(Event& ev, EventTmp& ev_tmp, MkBuilder& builder)
{
  EventOfCombCandidates &event_of_comb_cands = ev_tmp.m_event_of_comb_cands;
  event_of_comb_cands.Reset();

  builder.begin_event(&ev, &ev_tmp, __func__);

  if   (Config::findSeeds) {builder.find_seeds();}
  else                     {builder.map_seed_hits();} // all other simulated seeds need to have hit indices line up in LOH for seed fit

  builder.fit_seeds();

  builder.find_tracks_load_seeds();

#ifdef USE_VTUNE_PAUSE
  __itt_resume();
#endif

  double time = dtime();

  builder.FindTracksStandard();

  time = dtime() - time;

#ifdef USE_VTUNE_PAUSE
  __itt_pause();
#endif
  
  if (!Config::normal_val) {
    if (!Config::silent) builder.quality_output_COMB();
  } else {builder.root_val_COMB();}

  builder.end_event();

  return time;
}

//==============================================================================
// runBuildTestPlex Combinatorial: CloneEngine TBB
//==============================================================================

double runBuildingTestPlexCloneEngine(Event& ev, EventTmp& ev_tmp, MkBuilder& builder)
{
  EventOfCombCandidates &event_of_comb_cands = ev_tmp.m_event_of_comb_cands;
  event_of_comb_cands.Reset();

  builder.begin_event(&ev, &ev_tmp, __func__);

  if   (Config::findSeeds) {builder.find_seeds();}
  else                     {builder.map_seed_hits();} // all other simulated seeds need to have hit indices line up in LOH for seed fit

  builder.fit_seeds();

  builder.find_tracks_load_seeds();

#ifdef USE_VTUNE_PAUSE
  __itt_resume();
#endif
  double time = dtime();

  builder.FindTracksCloneEngine();

  time = dtime() - time;

#ifdef USE_VTUNE_PAUSE
  __itt_pause();
#endif

  if (!Config::normal_val) {
    if (!Config::silent) builder.quality_output_COMB();
  } else {builder.root_val_COMB();}

  builder.end_event();

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

  if   (Config::findSeeds) {builder.find_seeds();}
  else                     {builder.map_seed_hits();} // all other simulated seeds need to have hit indices line up in LOH for seed fit

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

  if (!Config::normal_val) {
    if (!Config::silent) builder.quality_output_COMB();
  } else {builder.root_val_COMB();}

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

    if   (Config::findSeeds) {builder.find_seeds();}
    else                     {builder.map_seed_hits();} // all other simulated seeds need to have hit indices line up in LOH for seed fit

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
    if (!Config::normal_val) {
      if (!Config::silent) builder.quality_output_BH(event_of_cands);
    } else {
      builder.root_val_BH(event_of_cands);
    }

    builder.end_event();
  }
  
  return time;
}
#endif
