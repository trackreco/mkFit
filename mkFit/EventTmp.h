#ifndef _eventtmp_h_

#include "HitStructures.h"
#include "MkFitter.h"
#include "CandCloner.h"

class MkFitter;
class CandCloner;

class EventTmp
{
public:
  EventOfCombCandidates m_event_of_comb_cands;

  // XXXX Should also add EventOfCandidates?
  // XXXX Should we go for delayed initialization/resizing then?

  std::vector<CandCloner*> m_cand_cloners;

  EventTmp();
  ~EventTmp();

  void AssureCandClonersExist(int n_thr);
  void DeleteCandCloners();
};

#endif
