#ifndef MkBuilder_h
#define MkBuilder_h

#include <vector>

//------------------------------------------------------------------------------

#include "MkFitter.h"
#include "CandCloner.h"
#include "SteeringParams.h"

#include <functional>
#include <mutex>

#include "Pool.h"

class TrackerInfo;
class LayerInfo;

//------------------------------------------------------------------------------

struct ExecutionContext
{
  Pool<CandCloner> m_cloners;
  Pool<MkFitter>   m_fitters;

  void populate(int n_thr)
  {
    m_cloners.populate(n_thr - m_cloners.size());
    m_fitters.populate(n_thr - m_fitters.size());
  }
};

extern ExecutionContext g_exe_ctx;

//==============================================================================
// The usual
//==============================================================================

class Event;
class EventTmp;

class MkBuilder
{
protected:
  void fit_one_seed_set(TrackVec& simtracks, int itrack, int end, MkFitter *mkfp,
                        const bool is_brl[], const SteeringParams &st_par);

  Event         *m_event;
  EventTmp      *m_event_tmp;
  EventOfHits    m_event_of_hits;

  int m_cnt=0, m_cnt1=0, m_cnt2=0, m_cnt_8=0, m_cnt1_8=0, m_cnt2_8=0, m_cnt_nomc=0;

  SteeringParams   m_steering_params[5];
  std::vector<int> m_brl_ecp_regions;

public:
  typedef std::vector<std::pair<int,int>> CandIdx_t;

  MkBuilder();
  ~MkBuilder();

  // --------

  static MkBuilder* make_builder();

  virtual void begin_event(Event* ev, EventTmp* ev_tmp, const char* build_type);

  int find_seeds();
  virtual void fit_seeds();

  void end_event();
  
  // --------

  void map_seed_hits(); // m_event->layerHits_ -> m_event_of_hits.m_layers_of_hits (seeds only)
  void remap_seed_hits(); // m_event_of_hits.m_layers_of_hits -> m_event->layerHits_ (seeds only)
  void remap_cand_hits(); // m_event_of_hits.m_layers_of_hits -> m_event->layerHits_ (cands only)
  void align_simtracks(); // simtrack labels get screwed up in endcap tests

  void quality_output_BH();
  void quality_output_COMB();
  void quality_reset();
  void quality_process(Track& tkcand);
  void quality_print();

  void quality_store_tracks_COMB();

  void root_val_BH();
  void root_val_COMB();
  void init_track_extras();

  // --------

  void find_tracks_load_seeds_BH(); // for FindTracksBestHit
  void find_tracks_load_seeds();
  void find_tracks_in_layers(EtaRegionOfCombCandidates &comb_candidates,
                             CandCloner &cloner, MkFitter *mkfp,
                             int start_seed, int end_seed, int region);

  // --------

  virtual void FindTracksBestHit();
  virtual void FindTracksStandard();
  virtual void FindTracksCloneEngine();
#ifdef USE_CUDA
  const Event* get_event() const { return m_event; }
  const EventOfHits& get_event_of_hits() const { return m_event_of_hits; }
#endif
};

#endif
