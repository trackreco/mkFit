#ifndef _buildtest_mplex_
#define _buildtest_mplex_

#include "Event.h"
#include "EventTmp.h"
#include "Track.h"

double runBuildingTestPlexBestHit(Event& ev);

double runBuildingTestPlexStandard(Event& ev, EventTmp& evtmp);

double runBuildingTestPlexCloneEngine(Event& ev, EventTmp& evtmp);

#if USE_CUDA
double runAllBuildingTestPlexBestHitGPU(std::vector<Event> &events);
#endif

#endif
