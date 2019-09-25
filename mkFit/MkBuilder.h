#ifndef MkBuilder_h
#define MkBuilder_h

#include <vector>
#include <map>

//------------------------------------------------------------------------------

#include "CandCloner.h"
#include "MkFitter.h"
#include "MkFinder.h"
#include "MkFinderFV.h"
#include "SteeringParams.h"

#include <functional>
#include <mutex>

#include "align_alloc.h"

#include "Pool.h"
//#define DEBUG
#include "Debug.h"

namespace mkfit {

class TrackerInfo;
class LayerInfo;

//------------------------------------------------------------------------------

using MkFinderFvVec = std::vector<MkFinderFv, aligned_allocator<MkFinderFv, 64>>;

struct ExecutionContext
{
  ExecutionContext() {}
  ~ExecutionContext() {
    dprint("MkFinderFvVec count " << m_finderv.unsafe_size());
  }

  Pool<CandCloner> m_cloners;
  Pool<MkFitter>   m_fitters;
  Pool<MkFinder>   m_finders;
  tbb::concurrent_queue<MkFinderFvVec> m_finderv;

  void populate(int n_thr)
  {
    m_cloners.populate(n_thr - m_cloners.size());
    m_fitters.populate(n_thr - m_fitters.size());
    m_finders.populate(n_thr - m_finders.size());
  }

  void populate_finderv(int n_thr, int n_seedsPerThread)
  {
    int count = MkFinderFv::nMplx(n_seedsPerThread);
    auto sz = m_finderv.unsafe_size();
    if (sz < static_cast<size_t>(n_thr)) {
      dprint("Allocating " << n_thr - sz << " MkFinderFvVec of length " << count);
      for (size_t i = 0; i < n_thr - sz; ++i) {
        MkFinderFvVec fv(count);
        m_finderv.push(std::move(fv));
      }
    }
  }

  MkFinderFvVec getFV(int sz)
  {
    size_t count = (sz + MkFinderFv::Seeds - 1)/MkFinderFv::Seeds; // adjust for seeds per finder
    MkFinderFvVec retfv;
    m_finderv.try_pop(retfv);
    if (retfv.size() < count) {
      dprint("Resizing from " << retfv.size() << " to " << count);
      retfv.resize(count);
    }
    return retfv;
  }

  void pushFV(MkFinderFvVec fv)
  {
    m_finderv.push(std::move(fv));
  }
};

extern ExecutionContext g_exe_ctx;

//==============================================================================
// The usual
//==============================================================================

class Event;

class MkBuilder
{
protected:
  void fit_one_seed_set(TrackVec& simtracks, int itrack, int end, MkFitter *mkfttr,
                        const bool is_brl[]);

  Event                 *m_event;
  EventOfHits            m_event_of_hits;
  EventOfCombCandidates  m_event_of_comb_cands;

  int m_cnt=0, m_cnt1=0, m_cnt2=0, m_cnt_8=0, m_cnt1_8=0, m_cnt2_8=0, m_cnt_nomc=0;

  FindingFoos      m_fndfoos_brl, m_fndfoos_ec;
  SteeringParams   m_steering_params[5];
  std::vector<int> m_regions;

public:
  typedef std::vector<std::pair<int,int>> CandIdx_t;

  MkBuilder();
  ~MkBuilder();

  // --------

  static MkBuilder* make_builder();
  static void populate(bool populatefv = false)
  {
    g_exe_ctx.populate(Config::numThreadsFinder);
    if (populatefv) {
      g_exe_ctx.populate_finderv(Config::numThreadsFinder, Config::numSeedsPerTask);
    }
  }

  int total_cands() const { 
    int res = 0; 
    for (auto const& icomb: m_event_of_comb_cands.m_candidates) res += icomb.size();
    return res;
  }

  std::pair<int,int> max_hits_layer() const {
    int maxN = 0;
    int maxL = 0;
    for (auto const& lh : m_event_of_hits.m_layers_of_hits){
      auto lsize = static_cast<int>(lh.m_hit_phis.size());
      if (lsize > maxN){
        maxN = lsize;
        maxL = lh.layer_id();
      }
    }
    return {maxN, maxL};
  }

  void begin_event(Event* ev, const char* build_type);
  void end_event();

  void create_seeds_from_sim_tracks();
  void import_seeds();
  void find_seeds();
  void assign_seedtype_forranking();
  void fit_seeds();

  // --------

  void map_track_hits  (TrackVec & tracks); // m_event->layerHits_ -> m_event_of_hits.m_layers_of_hits
  void remap_track_hits(TrackVec & tracks); // m_event_of_hits.m_layers_of_hits -> m_event->layerHits_

  void quality_val();
  void quality_reset();
  void quality_process(Track& tkcand, const int itrack, std::map<int,int> & cmsswLabelToPos);
  void quality_print();
  void track_print(Track &t, const char* pref);

  void quality_store_tracks(TrackVec & tracks);

  void root_val_dumb_cmssw();
  void root_val();
  void cmssw_export();
  void prep_recotracks(); 
  void prep_simtracks();
  void prep_cmsswtracks();
  void prep_reftracks(TrackVec& tracks, TrackExtraVec& extras, const bool realigntracks); 
  void prep_tracks(TrackVec& tracks, TrackExtraVec& extras, const bool realigntracks); // sort hits by layer, init track extras, align track labels if true
  void score_tracks(TrackVec& tracks); // if track score not already assigned

  void find_duplicates(TrackVec& tracks);
  void remove_duplicates(TrackVec& tracks);
  void handle_duplicates();

  // --------

  void find_tracks_load_seeds_BH(); // for FindTracksBestHit
  void find_tracks_load_seeds();

  int  find_tracks_unroll_candidates(std::vector<std::pair<int,int>> & seed_cand_vec,
                                     int start_seed, int end_seed,
                                     int prev_layer, bool pickup_only);

  void find_tracks_handle_missed_layers(MkFinder *mkfndr, const LayerInfo &layer_info,
                                        std::vector<std::vector<TrackCand>> &tmp_cands,
                                        const std::vector<std::pair<int,int>> &seed_cand_idx,
                                        const int region, const int start_seed,
                                        const int itrack, const int end);

  void find_tracks_in_layers(CandCloner &cloner, MkFinder *mkfndr,
                             const int start_seed, const int end_seed, const int region);
  void find_tracks_in_layersFV(int start_seed, int end_seed, int region);

  // --------

  void PrepareSeeds();

  void FindTracksBestHit();
  void FindTracksStandard();
  void FindTracksCloneEngine();
  void FindTracksFV();

  void BackwardFitBH();
  void fit_cands_BH(MkFinder *mkfndr, int start_cand, int end_cand, int region);

  void BackwardFit();
  void fit_cands(MkFinder *mkfndr, int start_cand, int end_cand, int region);

#ifdef USE_CUDA
  const Event* get_event() const { return m_event; }
  const EventOfHits& get_event_of_hits() const { return m_event_of_hits; }
#endif
};

} // end namespace mkfit
#endif
