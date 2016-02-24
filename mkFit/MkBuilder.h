#ifndef MkBuilder_h
#define MkBuilder_h

#include <vector>

#include "HitStructures.h"


//------------------------------------------------------------------------------

#include "MkFitter.h"
#include "CandCloner.h"

#include <stack>
#include <functional>
#include <mutex>

template <typename TT>
struct Pool
{
  typedef std::function<TT*()>     CFoo_t;
  typedef std::function<void(TT*)> DFoo_t;

  CFoo_t m_create_foo  = []()     { return new (_mm_malloc(sizeof(TT), 64)) TT; };
  DFoo_t m_destroy_foo = [](TT* x){ _mm_free(x); };

  std::stack<TT*> m_stack;
  std::mutex      m_moo;

  Pool() {}
  Pool(CFoo_t cf, DFoo_t df) : m_create_foo(cf), m_destroy_foo(df) {}

  ~Pool()
  {
    while ( ! m_stack.empty())
    {
      m_destroy_foo(m_stack.top());
      m_stack.pop();
    }
  }

  void SetCFoo(CFoo_t cf) { m_create_foo  = cf; }
  void SetDFoo(DFoo_t df) { m_destroy_foo = df; }

  TT* GetFromPool()
  {
    std::unique_lock<std::mutex> lk(m_moo);

    if (m_stack.empty()) return m_create_foo();

    TT *x = m_stack.top();
    m_stack.pop();
    return x;
  }

  void ReturnToPool(TT *x)
  {
    std::unique_lock<std::mutex> lk(m_moo);

    m_stack.push(x);
  }
};

struct ExecutionContext
{
  Pool<CandCloner> m_cloners;
  Pool<MkFitter>   m_fitters;

  ExecutionContext() {}
};


//==============================================================================
// The usual
//==============================================================================

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

  void FindTracksCloneEngineMT();
};

#endif
