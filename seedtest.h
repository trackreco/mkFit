#ifndef _seedtest_
#define _seedtest_

#include "Track.h"
#include "BinInfoUtils.h"
#include "Event.h"

void buildSeedsByMC(const TrackVec&, TrackVec&, Event&);
void buildSeedsByRoadTriplets(const std::vector<HitVec>&, const BinInfoMap&, TrackVec&, Event&);
void buildHitPairs(const std::vector<HitVec>&, const BinInfoLayerMap&, std::vector<HitVec>&);
void buildHitTriplets(const std::vector<HitVec>&, const BinInfoLayerMap&, const std::vector<HitVec>&, std::vector<HitVec>&);
void filterHitTripletsByRZChi2(const std::vector<HitVec> &, std::vector<HitVec>&);
void buildSeedsFromTriplets(const std::vector<HitVec>&, TrackVec&, Event&);

#endif
