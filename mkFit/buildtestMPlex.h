#ifndef _buildtest_mplex_
#define _buildtest_mplex_

#include "Event.h"
#include "EventTmp.h"
#include "Track.h"

double runBuildingTestPlexBestHit(Event& ev);

double runBuildingTestPlex(Event& ev, EventTmp& ev_tmp);

double runBuildingTestPlexCloneEngine(Event& ev, EventTmp& evtmp);

double runBuildingTestPlexTbb(Event& ev, EventTmp& evtmp);

#endif
