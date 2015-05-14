#ifndef _seedtest_
#define _seedtest_

#include "Track.h"
#include "BinInfoUtils.h"

void buildSeedsByMC(const TrackVec&, TrackVec&);
void buildSeedsByRoadTriplets(const std::vector<HitVec>&, const BinInfoMap&, TrackVec&);
void buildHitPairs(const std::vector<HitVec>&, const BinInfoLayerMap&, std::vector<HitVec>&);
void buildHitTriplets(const std::vector<HitVec>&, const BinInfoLayerMap&, std::vector<HitVec>&);
void filterHitTripletsByRZChi2(std::vector<HitVec> &);
void buildSeedsFromTriplets(const std::vector<HitVec>&, TrackVec&);

#endif
