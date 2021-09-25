#ifndef MkBuilder_h
#define MkBuilder_h

#include <vector>
#include <map>

//------------------------------------------------------------------------------

#include "CandCloner.h"
#include "MkFitter.h"
#include "MkFinder.h"
#include "IterationConfig.h"

#include <functional>
#include <mutex>

#include "align_alloc.h"

#include "Pool.h"
//#define DEBUG
#include "Debug.h"

namespace mkfit {

class TrackerInfo;
class LayerInfo;

class Event;

//==============================================================================
// Execution context -- Pools of helper objects
//==============================================================================

struct ExecutionContext
{
  ExecutionContext() = default;
  ~ExecutionContext() = default;

  Pool<CandCloner> m_cloners;
  Pool<MkFitter>   m_fitters;
  Pool<MkFinder>   m_finders;

  void populate(int n_thr)
  {
    m_cloners.populate(n_thr - m_cloners.size());
    m_fitters.populate(n_thr - m_fitters.size());
    m_finders.populate(n_thr - m_finders.size());
  }
};

extern ExecutionContext g_exe_ctx;

//==============================================================================
// MkJob
//==============================================================================

class MkJob
{
public:
  const TrackerInfo          &m_trk_info;
  // Config &config; // If we want to get rid of namespace / global config
  const IterationConfig      &m_iter_config;
  const EventOfHits          &m_event_of_hits;

  const IterationMaskIfcBase *m_iter_mask_ifc = nullptr;

        int  num_regions()   const { return m_iter_config.m_n_regions; }
  const auto regions_begin() const { return m_iter_config.m_region_order.begin(); }
  const auto regions_end()   const { return m_iter_config.m_region_order.end(); }

  const auto& steering_params(int i) { return m_iter_config.m_steering_params[i]; }

  const auto& params()     const { return m_iter_config.m_params; }
  const auto& params_bks() const { return m_iter_config.m_backward_params; }

  int max_max_cands() const
  { return std::max(params().maxCandsPerSeed, params_bks().maxCandsPerSeed); }

  const std::vector<bool>* get_mask_for_layer(int layer)
  {
    return m_iter_mask_ifc ? m_iter_mask_ifc->get_mask_for_layer(layer) : nullptr;
  }
};


//==============================================================================
// MkBuilder
//==============================================================================

class MkBuilder
{
protected:
  void fit_one_seed_set(TrackVec& simtracks, int itrack, int end, MkFitter *mkfttr,
                        const bool is_brl[]);

  MkJob     *m_job  = nullptr;

  // MIMI -- To be removed. Used by seed processing / validation that has yet to be moved.
  Event     *m_event = nullptr;

  // State for BestHit
  TrackVec                        m_tracks;

  // State for Std / CloneEngine
  EventOfCombCandidates           m_event_of_comb_cands;

  int m_cnt=0, m_cnt1=0, m_cnt2=0, m_cnt_8=0, m_cnt1_8=0, m_cnt2_8=0, m_cnt_nomc=0;

  // Per-region seed information
  IntVec           m_seedEtaSeparators;
  IntVec           m_seedMinLastLayer;
  IntVec           m_seedMaxLastLayer;

  std::atomic<int> m_nan_n_silly_per_layer_count;

public:
  using insert_seed_foo       = void(const Track &, int);
  using filter_track_cand_foo = bool(const TrackCand &);

  typedef std::vector<std::pair<int,int>> CandIdx_t;

  MkBuilder();
  ~MkBuilder();

  // --------

  static MkBuilder* make_builder();
  static void populate()
  {
    g_exe_ctx.populate(Config::numThreadsFinder);
  }

  int total_cands() const
  {
    int res = 0;
    for (auto const& icomb: m_event_of_comb_cands.m_candidates) res += icomb.size();
    return res;
  }

  std::pair<int,int> max_hits_layer(const EventOfHits &eoh) const
  {
    int maxN = 0;
    int maxL = 0;
    for (auto const& lh : eoh.m_layers_of_hits)
    {
      auto lsize = static_cast<int>(lh.m_hit_phis.size());
      if (lsize > maxN){
        maxN = lsize;
        maxL = lh.layer_id();
      }
    }
    return {maxN, maxL};
  }

  void begin_event(MkJob *job, Event *ev, const char *build_type);
  void end_event();
  void import_seeds(const TrackVec &in_seeds, std::function<insert_seed_foo> insert_seed);

  int  filter_comb_cands(std::function<filter_track_cand_foo> filter);

  // XXX filter for rearranging cands that will / will not do backward search.

  void select_best_comb_cands(bool clear_m_tracks=false);
  void export_best_comb_cands(TrackVec &out_vec);
  void export_tracks(TrackVec &out_vec);

  void CompactifyHitStorageForBestCand(bool remove_seed_hits, int backward_fit_min_hits)
  { m_event_of_comb_cands.CompactifyHitStorageForBestCand(remove_seed_hits, backward_fit_min_hits); }

  void BeginBkwSearch() { m_event_of_comb_cands.BeginBkwSearch(); }
  void EndBkwSearch()   { m_event_of_comb_cands.EndBkwSearch(); }

  // MIMI hack to export tracks for BH
  const TrackVec& ref_tracks() const { return m_tracks; }
        TrackVec& ref_tracks_nc()    { return m_tracks; }

  // void create_seeds_from_sim_tracks();
  // void find_seeds();
  // void fit_seeds();

  // --------

  void quality_val();
  void quality_reset();
  void quality_process(Track& tkcand, const int itrack, std::map<int,int> & cmsswLabelToPos);
  void quality_print();
  void track_print(const Track &t, const char* pref);

  void quality_store_tracks(TrackVec & tracks);

  void root_val_dumb_cmssw();
  void root_val();

  void prep_recotracks();
  void prep_simtracks();
  void prep_cmsswtracks();
  void prep_reftracks(TrackVec& tracks, TrackExtraVec& extras, const bool realigntracks);
  void prep_tracks(TrackVec& tracks, TrackExtraVec& extras, const bool realigntracks); // sort hits by layer, init track extras, align track labels if true
  void score_tracks(TrackVec& tracks); // if track score not already assigned

  // --------

  void find_tracks_load_seeds_BH(const TrackVec &in_seeds); // for FindTracksBestHit
  void find_tracks_load_seeds   (const TrackVec &in_seeds);

  int  find_tracks_unroll_candidates(std::vector<std::pair<int,int>> & seed_cand_vec,
                                     int start_seed, int end_seed,
                                     int layer, int prev_layer, bool pickup_only,
                                     SteeringParams::IterationType_e iteration_dir);

  void find_tracks_handle_missed_layers(MkFinder *mkfndr, const LayerInfo &layer_info,
                                        std::vector<std::vector<TrackCand>> &tmp_cands,
                                        const std::vector<std::pair<int,int>> &seed_cand_idx,
                                        const int region, const int start_seed,
                                        const int itrack, const int end);

  void find_tracks_in_layers(CandCloner &cloner, MkFinder *mkfndr,
                             SteeringParams::IterationType_e iteration_dir,
                             const int start_seed, const int end_seed, const int region);

  // --------
  static void seed_post_cleaning(TrackVec &tv, const bool fix_silly_seeds, const bool remove_silly_seeds);

  void PrepareSeeds();

  void FindTracksBestHit(SteeringParams::IterationType_e iteration_dir=SteeringParams::IT_FwdSearch);
  void FindTracksStandard(SteeringParams::IterationType_e iteration_dir=SteeringParams::IT_FwdSearch);
  void FindTracksCloneEngine(SteeringParams::IterationType_e iteration_dir=SteeringParams::IT_FwdSearch);

  void BackwardFitBH();
  void fit_cands_BH(MkFinder *mkfndr, int start_cand, int end_cand, int region);

  void BackwardFit();
  void fit_cands(MkFinder *mkfndr, int start_cand, int end_cand, int region);
};

} // end namespace mkfit
#endif
