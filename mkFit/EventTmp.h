#ifndef _eventtmp_h_

#include "HitStructures.h"
#include "MkFitter.h"

class MkFitter;

class EventTmp
{
public:
  EventOfCombCandidates m_event_of_comb_cands;

  // XXXX Should also add EventOfCandidates?
  // XXXX Should we go for delayed initialization/resizing then?

  EventTmp() {}
  ~EventTmp() {}
};

#endif
