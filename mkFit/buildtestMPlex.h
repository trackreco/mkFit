#ifndef _buildtest_mplex_
#define _buildtest_mplex_

#include "Event.h"
#include "Track.h"
#ifdef USE_CUDA
#include "BuilderCU.h"
#endif

class MkBuilder;

double runBuildingTestPlexBestHit(Event& ev, MkBuilder& builder);
double runBuildingTestPlexStandard(Event& ev, MkBuilder& builder);
double runBuildingTestPlexCloneEngine(Event& ev, MkBuilder& builder);
double runBuildingTestPlexFV(Event& ev, MkBuilder& builder);

#if USE_CUDA
double runBuildingTestPlexBestHitGPU(Event& ev, MkBuilder& builder,
                                     BuilderCU& builder_cu);
double runBuildingTestPlexCloneEngineGPU(Event& ev, EventTmp& ev_tmp, MkBuilder& builder,
                                         BuilderCU& builder_cu, bool seed_based=false);
double runAllBuildingTestPlexBestHitGPU(std::vector<Event> &events);
#endif

#endif
