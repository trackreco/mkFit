#include "buildtestMPlex.h"
#include "Matrix.h"
#include "MkBuilder.h"
#include "MkStdSeqs.h"

#include "tbb/parallel_for.h"

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
    // MIMI -- nan_n_silly_per_layer_count is in MkBuilder, could be in MkJob.
    // if (Config::nan_n_silly_check_cands_every_layer)
    // {
    //   int sc = (int) ev.nan_n_silly_per_layer_count_;
    //   if (sc > 0)
    //     printf("Nan'n'Silly: Number of silly candidates over all layers = %d\n", sc);
    // }
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

void runBuildingTestPlexDumbCMSSW(Event& ev, const EventOfHits &eoh, MkBuilder& builder)
{
  MkJob job( { Config::TrkInfo, Config::ItrInfo[0], eoh } );

  builder.begin_event(&job, &ev, __func__);

  if (Config::sim_val_for_cmssw) {
    builder.root_val_dumb_cmssw();
  }

  builder.end_event();
}

//==============================================================================
// runBuildTestPlexBestHit
//==============================================================================

double runBuildingTestPlexBestHit(Event& ev, const EventOfHits &eoh, MkBuilder& builder)
{
  MkJob job( { Config::TrkInfo, Config::ItrInfo[0], eoh } );

  builder.begin_event(&job, &ev, __func__);

  builder.PrepareSeeds();

  // EventOfCandidates event_of_cands;
  builder.find_tracks_load_seeds_BH(ev.seedTracks_);

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

  // Hack, get the tracks out.
  ev.candidateTracks_ = builder.ref_tracks();

  // For best hit, the candidateTracks_ vector is the direct input to the backward fit so only need to do find_duplicates once
  if (Config::quality_val || Config::sim_val || Config::cmssw_val || Config::cmssw_export)
  {
    //Mark tracks as duplicates; if within CMSSW, remove duplicate tracks before backward fit
    if(Config::removeDuplicates)
    {
      StdSeq::find_duplicates(ev.candidateTracks_);
      if(Config::cmssw_export) StdSeq::remove_duplicates(ev.candidateTracks_);
    }
  }

  // now do backwards fit... do we want to time this section?
  if (Config::backwardFit)
  {
    builder.BackwardFitBH();
    ev.fitTracks_ = builder.ref_tracks();
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

double runBuildingTestPlexStandard(Event& ev, const EventOfHits &eoh, MkBuilder& builder)
{
  MkJob job( { Config::TrkInfo, Config::ItrInfo[0], eoh } );

  builder.begin_event(&job, &ev, __func__);

  builder.PrepareSeeds();

  builder.find_tracks_load_seeds(ev.seedTracks_);

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
    // Using the TrackVec version until we home in on THE backward fit etc.
    // builder.BackwardFit();
    builder.select_best_comb_cands();
    builder.BackwardFitBH();
    ev.fitTracks_ = builder.ref_tracks();

    check_nan_n_silly_bkfit(ev);
  }

  StdSeq::handle_duplicates(&ev);

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

double runBuildingTestPlexCloneEngine(Event& ev, const EventOfHits &eoh, MkBuilder& builder)
{
  MkJob job( { Config::TrkInfo, Config::ItrInfo[0], eoh } );

  builder.begin_event(&job, &ev, __func__);

  builder.PrepareSeeds();

  builder.find_tracks_load_seeds(ev.seedTracks_);

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

  // first store candidate tracks - needed for BH backward fit and root_validation
  builder.quality_store_tracks(ev.candidateTracks_);

  // now do backwards fit... do we want to time this section?
  if (Config::backwardFit)
  {
    // a) TrackVec version:
    builder.select_best_comb_cands();
    builder.BackwardFitBH();
    ev.fitTracks_ = builder.ref_tracks();

    // b) Version that runs on CombCand / TrackCand
    // builder.BackwardFit();
    // builder.quality_store_tracks(ev.fitTracks_);

    check_nan_n_silly_bkfit(ev);
  }
  
  StdSeq::handle_duplicates(&ev);
  
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
// runBtbCe_MultiIter
//
// Prototype for running multiple iterations, sequentially, using the same builder.
// For cmmsw seeds
//
// There is, in general, a mess in how tracks are processed, marked, or copied out
// in various validation scenarios and export flags.
//
// In particular, MkBuilder::PrepareSeeds does a lot of things to whole / complete
// event,seedTracks_ -- probably this would need to be split into common / and
// per-iteration part.
// - hit remapping back happens in MkBuilder::quality_val() and root_validation().
// - MkBuilder::prep_*** functions also mostly do not belong there (prep_sim is called from
//   PrepareSeeds() for cmssw seeds).
//
// StdSeq::handle_duplicates() also goes over all Event tracks ... should it be per iteration?
//
// At this point we need to think about what should happen to Event before all the iterations and
// after all the iterations ... from the Validation perspective.
// And if we care about doing too muich work for seeds that will never get processed.
//==============================================================================

double runBtbCe_MultiIter(Event& ev, const EventOfHits &eoh, MkBuilder& builder, int n)
{
  double ttime = 0;
  if (n<=0) return ttime; //at least one iter by default

  const bool validation_on = (Config::sim_val || Config::quality_val);
  
  TrackVec seeds_used;
  TrackVec seeds1;
   
  unsigned int algorithms[]={ 4,22,23,5,24,7,8,9,10 }; //9 iterations

  if (validation_on) 
  {
    for (auto &s : ev.seedTracks_)
    {
      //keep seeds form the first n iterations for processing
      if ( std::find(algorithms, algorithms+n, s.algoint())!=algorithms+n  ) seeds1.push_back(s);
    }
    ev.seedTracks_.swap(seeds1);//necessary for the validation - PrepareSeeds
    ev.relabel_bad_seedtracks();//necessary for the validation - PrepareSeeds
  }
  
  IterationMaskIfc mask_ifc;

  for (int it = 0; it <= n-1; ++it)
  {
    // MIMI - to disable hit-masks, pass nullptr in place of &mask_ifc to job
    // and optionally comment out ev.fill_hitmask_bool_vectors() call.

    ev.fill_hitmask_bool_vectors(Config::ItrInfo[it].m_track_algorithm, mask_ifc.m_mask_vector);

    MkJob job( { Config::TrkInfo, Config::ItrInfo[it], eoh, &mask_ifc } );

    builder.begin_event(&job, &ev, __func__);

    // Some of what happens here, should really happen somewhere else, if it needs to.
    // builder.PrepareSeeds();
    // Specific cleaning / mapping done below for extracted seeds.

    TrackVec seeds;
    { // We could partition seeds once, store beg, end for each iteration in a map or vector.
      int nc = 0;
      for (auto &s : ev.seedTracks_)
      {
        if (s.algoint() == job.m_iter_config.m_track_algorithm) {
          seeds.push_back(s);
          ++nc;
        }
        else if (nc > 0) break;
      }
    }

    // MIMI -- using Iter0 function / tuning for all iterations.
    ev.clean_cms_seedtracks_iter(&seeds, Config::ItrInfo[it]);

    if(seeds.size()<=0)
      continue;
    
    builder.seed_post_cleaning(seeds, true, true);

    for (auto &s : seeds) assignSeedTypeForRanking(s);

    builder.find_tracks_load_seeds(seeds);

    double time = dtime();

    builder.FindTracksCloneEngine();

    ttime += dtime() - time;

    if (validation_on)  seeds_used.insert(seeds_used.end(), seeds.begin(), seeds.end());//cleaned seeds need to be stored somehow
    
    // first store candidate tracks - needed for BH backward fit and root_validation
    // XXXX to builder m_tracks ... or do we do this for validation anyway ?    
    builder.export_best_comb_cands(ev.candidateTracks_);

    // now do backwards fit... do we want to time this section?
    if (Config::backwardFit)
    {
      // a) TrackVec version:
      builder.select_best_comb_cands();
      builder.BackwardFitBH();
      builder.export_tracks(ev.fitTracks_);

      // b) Version that runs on CombCand / TrackCand
      // builder.BackwardFit();
      // builder.quality_store_tracks(ev.fitTracks_);

    }
    
    builder.end_event();
  }

  // MIMI - Fake back event pointer for final processing (that should be done elsewhere)
  MkJob job( { Config::TrkInfo, Config::ItrInfo[0], eoh } );
  builder.begin_event(&job, &ev, __func__);
  
  //I have the seeds -> continue the PrepareSeeds()
  if (validation_on)
  {
      builder.prep_simtracks();
      //swap for the cleaned seeds
      ev.seedTracks_.swap(seeds_used);
  }    
  //PrepareSeeds() done
  
  check_nan_n_silly_candiates(ev);

  if (Config::backwardFit) check_nan_n_silly_bkfit(ev);

  // ??? MIMI - How to do this? Assuming cleaning of all iterations together.
  StdSeq::handle_duplicates(&ev);

  // ??? MIMI - And same here ... most validation runs on Event, I hope.
  // validation section
  if        (Config::quality_val) {
    builder.quality_val();
  } else if (Config::sim_val || Config::cmssw_val) {
    builder.root_val();
  } else if (Config::cmssw_export) {
    builder.cmssw_export();
  }

  // ev.print_tracks(ev.candidateTracks_, true);

  // MIMI Unfake.
  builder.end_event();

  return ttime;
}


//==============================================================================
// run_OneIteration
//
// One-stop function for running track building from CMSSW.
//==============================================================================

struct IterationMaskIfcCmssw : public IterationMaskIfcBase
{
  const TrackerInfo                           &m_trk_info;
  const std::vector<const std::vector<bool>*> &m_mask_vector;

  IterationMaskIfcCmssw(const TrackerInfo &ti, const std::vector<const std::vector<bool>*> &maskvec) :
    m_trk_info(ti), m_mask_vector(maskvec) {}

  const std::vector<bool>* get_mask_for_layer(int layer) const
  {
     return m_trk_info.m_layers[layer].is_pix_lyr() ? m_mask_vector[0] : m_mask_vector[1];
  }
};

void run_OneIteration(const TrackerInfo& trackerInfo, const IterationConfig &itconf, const EventOfHits &eoh,
                      const std::vector<const std::vector<bool>*>& hit_masks,
                      MkBuilder& builder, TrackVec &seeds, TrackVec &out_tracks,
                      bool do_seed_clean, bool do_backward_fit, bool do_remove_duplicates)
{
  IterationMaskIfcCmssw it_mask_ifc(trackerInfo, hit_masks);

  MkJob job( { trackerInfo, itconf, eoh, &it_mask_ifc } );

  builder.begin_event(&job, nullptr, __func__);

  if (do_seed_clean)
  {
  // XXXX MIMI -- Leonardo working on iter-dependent seed cleaning, will be out of Event!
    Event dummy_ev(-1);
    dummy_ev.clean_cms_seedtracks(&seeds);
  }

  // Check nans in seeds -- this should not be needed when Slava fixes
  // the track parameter coordinate transformation.
  builder.seed_post_cleaning(seeds, true, true);

  for (auto &s : seeds) assignSeedTypeForRanking(s);

  builder.find_tracks_load_seeds(seeds);

  builder.FindTracksCloneEngine();

  if ( ! do_backward_fit)
  {
    builder.export_best_comb_cands(out_tracks);
  }
  else
  {
    // a) BackwardFitBH works on m_tracks
    builder.select_best_comb_cands();
    builder.BackwardFitBH();
    builder.export_tracks(out_tracks);

    // b) Version that runs on CombCand / TrackCand
    // builder.BackwardFit();
    // builder.export_best_comb_cands(out_tracks);
  }

  if (do_remove_duplicates)
  {
    StdSeq::find_duplicates(out_tracks);
    StdSeq::remove_duplicates(out_tracks);
  }

  builder.end_event();
}

} // end namespace mkfit
