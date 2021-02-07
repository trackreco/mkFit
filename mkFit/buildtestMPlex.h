#ifndef _buildtest_mplex_
#define _buildtest_mplex_

#include "Event.h"
#include "Track.h"

#include "HitStructures.h"

namespace mkfit {

class MkBuilder;

void   runBuildingTestPlexDumbCMSSW  (Event& ev, EventOfHits &eoh, MkBuilder& builder);
double runBuildingTestPlexBestHit    (Event& ev, EventOfHits &eoh, MkBuilder& builder);
double runBuildingTestPlexStandard   (Event& ev, EventOfHits &eoh, MkBuilder& builder);
double runBuildingTestPlexCloneEngine(Event& ev, EventOfHits &eoh, MkBuilder& builder);

double runBtbCe_MultiIter(Event& ev, EventOfHits &eoh, MkBuilder& builder);

} // end namespace mkfit
#endif
