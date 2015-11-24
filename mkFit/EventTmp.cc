#include "EventTmp.h"

EventTmp::EventTmp() :
  m_event_of_comb_cands()
{}

EventTmp::~EventTmp()
{
  DeleteCandCloners();
}

//------------------------------------------------------------------------------

void EventTmp::AssureCandClonersExist(int n_thr)
{
  if (m_cand_cloners.size() == n_thr) return;

  DeleteCandCloners();

  m_cand_cloners.resize(n_thr);

  for (int i = 0; i < n_thr; ++i)
  {
    if (Config::clonerUseSingleThread)
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
  }
}

void EventTmp::DeleteCandCloners()
{
  int n_thr = m_cand_cloners.size();

  for (int i = 0; i < n_thr; ++i)
  {
    delete m_cand_cloners[i];
  }
}
