#ifndef MkBuilder_h
#define MkBuilder_h

#include <vector>

#include "HitStructures.h"


//------------------------------------------------------------------------------

#include "MkFitter.h"
#include "CandCloner.h"

#include "tbb/concurrent_queue.h"
#include <functional>
#include <mutex>

template <typename TT>
struct Pool
{
  typedef std::function<TT*()>     CFoo_t;
  typedef std::function<void(TT*)> DFoo_t;

  CFoo_t m_create_foo  = []()     { return new (_mm_malloc(sizeof(TT), 64)) TT; };
  DFoo_t m_destroy_foo = [](TT* x){ _mm_free(x); };

  tbb::concurrent_queue<TT*> m_stack;

  void populate()
  {
    for (int i = 0; i < Config::numThreadsFinder; ++i)
    {
      m_stack.push(m_create_foo());
    }
  }

  Pool() {}
  Pool(CFoo_t cf, DFoo_t df) : m_create_foo(cf), m_destroy_foo(df) {}

  ~Pool()
  {
    TT *x;
    while (m_stack.try_pop(x))
    {
      m_destroy_foo(x);
    }
  }

  void SetCFoo(CFoo_t cf) { m_create_foo  = cf; }
  void SetDFoo(DFoo_t df) { m_destroy_foo = df; }

  TT* GetFromPool()
  {
    TT *x;
    if (m_stack.try_pop(x)) {
      return x;
    } else {
      return m_create_foo();
    }
  }

  void ReturnToPool(TT *x)
  {
    m_stack.push(x);
  }
};

struct ExecutionContext
{
  Pool<CandCloner> m_cloners;
  Pool<MkFitter>   m_fitters;

  ExecutionContext()
  {
    m_cloners.populate();
    m_fitters.populate();
  }
};

extern ExecutionContext g_exe_ctx;

//==============================================================================
// The usual
//==============================================================================

class Event;
class EventTmp;

class MkFitter;

class MkBuilder
{
protected:
  void fit_one_seed_set(TrackVec& simtracks, int itrack, int end, MkFitter *mkfp);

  Event         *m_event;
  EventTmp      *m_event_tmp;
  EventOfHits    m_event_of_hits;

  std::vector<MkFitter*> m_mkfp_arr;

  int m_cnt=0, m_cnt1=0, m_cnt2=0, m_cnt_8=0, m_cnt1_8=0, m_cnt2_8=0, m_cnt_nomc=0;

public:
  typedef std::vector<std::pair<int,int>> CandIdx_t;

  MkBuilder();
  ~MkBuilder();

  // --------

  virtual void begin_event(Event* ev, EventTmp* ev_tmp, const char* build_type);

  virtual void fit_seeds();
  virtual void fit_seeds_tbb();

  void end_event();

  // --------

  void quality_reset();
  void quality_process(Track& tkcand);
  void quality_print();

  // --------

  // Common foos for FindTracks() / FindTracksCloneEngine() ???

  void find_tracks_load_seeds(EventOfCandidates& event_of_cands); // for FindTracksBestHit
  void find_tracks_load_seeds();
  void find_tracks_in_layers(EtaBinOfCombCandidates &eb_of_cc, CandCloner &cloner, MkFitter *mkfp,
                             int start_seed, int end_seed, int ebin);

  // --------

  virtual void FindTracksBestHit(EventOfCandidates& event_of_cands);
  virtual void FindTracks();
  virtual void FindTracksCloneEngine();
  virtual void FindTracksCloneEngineTbb();
};

#endif
