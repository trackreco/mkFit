#ifndef _seedtest_mplex_
#define _seedtest_mplex_

#include "Event.h"
#include "Track.h"
#include "HitStructures.h"

namespace mkfit {

void findSeedsByRoadSearch(TripletIdxConVec & seed_idcs, std::vector<LayerOfHits>& evt_lay_hits, int lay1_size, Event *& ev);

} // end namespace mkfit
#endif
