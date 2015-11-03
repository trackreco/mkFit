#include "EventTmp.h"

EventTmp::EventTmp() :
  m_event_of_comb_cands()
{
  // NOTE: MkFitter *MUST* be on heap, not on stack!
  // Standard operator new screws up alignment of ALL MPlex memebrs of MkFitter,
  // even if one adds attr(aligned(64)) thingy to every possible place.

  m_cand_cloners.resize(Config::g_num_threads);

  for (int i = 0; i < Config::g_num_threads; ++i)
  {
#ifdef TEST_CLONE_ENGINE

    if (Config::g_cloner_single_thread)
    {
#if defined(__MIC__)
      m_cand_cloners[i] = new CandCloner(1 + i*4, -1, false);
#else
      m_cand_cloners[i] = new CandCloner(i, -1, false);
#endif
    } else {
#if defined(__MIC__)
      // Same core
      //m_cand_cloners[i] = new CandCloner(1 + i*4, 1 + i*4 + 1, false);
      // Dedicated core
      m_cand_cloners[i] = new CandCloner(1 + 2*i*4, 1 + (2*i + 1)*4, false);
#else
      m_cand_cloners[i] = new CandCloner(2*i, 2*i + 1, false);
      // CandCloner cloner(8, 20); // Same core
      // CandCloner cloner(1, 2);  // Same socket, another core
      // CandCloner cloner(1, 7);  // Another socket
#endif
    }
#endif
  }
}

EventTmp::~EventTmp()
{
  for (int i = 0; i < Config::g_num_threads; ++i)
  {
#ifdef TEST_CLONE_ENGINE
    delete m_cand_cloners[i];
#endif
  }
}
