#ifndef _seedtest_mplex_
#define _seedtest_mplex_

#include "Event.h"
#include "Track.h"
#include "HitStructures.h"

void findSeedsByRoadSearch(TrackVec& evt_seed_tracks, std::vector<LayerOfHits>& evt_lay_hits, Event *& ev);

#endif
