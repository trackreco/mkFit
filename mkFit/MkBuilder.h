#ifndef MkBuilder_h
#define MkBuilder_h

#include <vector>

#include "HitStructures.h"
#include "FitterCU.h"

class Event;
class EventTmp;

class MkFitter;

class MkBuilder
{
private:
  Event         *m_event;
  EventTmp      *m_event_tmp;
  EventOfHits    m_event_of_hits;

  std::vector<MkFitter*> m_mkfp_arr;
#ifdef USE_CUDA
  std::vector<FitterCU<float>*> m_cuFitter_arr;
#endif

  std::vector<Track>     m_recseeds;

  int m_cnt=0, m_cnt1=0, m_cnt2=0, m_cnt_8=0, m_cnt1_8=0, m_cnt2_8=0, m_cnt_nomc=0;

public:
  MkBuilder();
  ~MkBuilder();

  // --------

  void begin_event(Event* ev, EventTmp* ev_tmp, const char* build_type);

  void fit_seeds();

  void end_event();

  // --------

  void quality_reset();
  void quality_process(Track& tkcand);
  void quality_print();

  // --------

  // Common foos for FindTracks() / FindTracksCloneEngine() ???

  void find_tracks_load_seeds();

  // --------

  void FindTracksBestHit(EventOfCandidates& event_of_cands);

  void FindTracks();

  void FindTracksCloneEngine();
};

#endif
