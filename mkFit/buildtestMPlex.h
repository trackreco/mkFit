#ifndef _buildtest_mplex_
#define _buildtest_mplex_

#include "Event.h"
#include "Track.h"

class MkBuilder;

double runBuildingTestPlexBestHit(Event& ev, MkBuilder& builder);
double runBuildingTestPlexStandard(Event& ev, MkBuilder& builder);
double runBuildingTestPlexCloneEngine(Event& ev, MkBuilder& builder);

#if USE_CUDA
double runAllBuildingTestPlexBestHitGPU(std::vector<Event> &events);
#endif

#endif
