#include "buildtestMPlex.h"

#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"

#include "MkBuilder.h"

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

  double time = dtime();

  builder.FindTracksBestHit();

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
  builder.quality_store_tracks(ev.candidateTracks_);

  // now do backwards fit... do we want to time this section?
  if (Config::backwardFit)
  {
    // QQQQ Using the TrackVec version until we home in on THE backward fit etc.
    // builder.BackwardFit();
    builder.BackwardFitBH();

    check_nan_n_silly_bkfit(ev);

    // QQQQ already done by BackwardFitBH()
    // if (Config::sim_val || Config::cmssw_val || Config::cmssw_export)
    // {
    //   builder.quality_store_tracks(ev.fitTracks_);
    // }
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
  builder.quality_store_tracks(ev.candidateTracks_);

  // now do backwards fit... do we want to time this section?
  if (Config::backwardFit)
  {
    // QQQQ Using the TrackVec version until we home in on THE backward fit etc.
    // builder.BackwardFit();
    builder.BackwardFitBH();

    check_nan_n_silly_bkfit(ev);

   // QQQQ already done by BackwardFitBH()
   // if (Config::sim_val || Config::cmssw_val || Config::cmssw_export)
   //  {
   //    builder.quality_store_tracks(ev.fitTracks_);
   //  }
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
} // end namespace mkfit
