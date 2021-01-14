#ifndef _buildtest_mplex_
#define _buildtest_mplex_

#include "Event.h"
#include "Track.h"

namespace mkfit {

class MkBuilder;

void   runBuildingTestPlexDumbCMSSW(Event& ev, MkBuilder& builder);
double runBuildingTestPlexBestHit(Event& ev, MkBuilder& builder);
double runBuildingTestPlexStandard(Event& ev, MkBuilder& builder);
double runBuildingTestPlexCloneEngine(Event& ev, MkBuilder& builder);
double runBuildingTestPlexFV(Event& ev, MkBuilder& builder);

} // end namespace mkfit
#endif
