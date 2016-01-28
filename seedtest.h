#ifndef _seedtest_
#define _seedtest_

#include "Track.h"
#include "BinInfoUtils.h"
#include "Event.h"
#include "Hit.h"

void buildSeedsByMC(const TrackVec&, TrackVec&, TrackExtraVec&, Event&);
void buildSeedsByRoadTriplets(TrackVec&, TrackExtraVec&, const std::vector<HitVec>&, const BinInfoMap&, Event&);
void buildHitPairs(const std::vector<HitVec>&, const BinInfoLayerMap&, std::vector<std::vector<int> >&);
void buildHitTriplets(const std::vector<HitVec>&, const BinInfoLayerMap&, const std::vector<std::vector<int> >&, std::vector<std::vector<int> >&);
void filterHitTripletsByRZChi2(const std::vector<HitVec>&, const std::vector<std::vector<int> >&, std::vector<std::vector<int> >&);
void buildSeedsFromTriplets(const std::vector<HitVec>&, const std::vector<std::vector<int> >&, TrackVec&, TrackExtraVec&, Event&);
int   getCharge(const Hit &, const Hit &, const Hit &);
void  rotateHitPos(float&, float&, const float);
float extrapolateToInterceptX(const float, const float, const float, const float, const float, const float);

#endif
