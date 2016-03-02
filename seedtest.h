#ifndef _seedtest_
#define _seedtest_

#include "Track.h"
#include "BinInfoUtils.h"
#include "Event.h"
#include "Hit.h"

void buildSeedsByMC(const TrackVec&, TrackVec&, TrackExtraVec&, Event&);
void buildSeedsByRZFirstRPhiSecond(TrackVec&, TrackExtraVec&, const std::vector<HitVec>&, const BinInfoMap&, Event&);
void filterHitTripletsByCircleParams(const std::vector<HitVec>&, const TripletIdxVec&, TripletIdxVec&);
void buildSeedsByRoadTriplets(TrackVec&, TrackExtraVec&, const std::vector<HitVec>&, const BinInfoMap&, Event&);
void buildHitPairs(const std::vector<HitVec>&, const BinInfoLayerMap&, PairIdxVec&);
void buildHitTripletsCurve(const std::vector<HitVec>&, const BinInfoLayerMap&, const PairIdxVec&, TripletIdxVec&);
void buildHitTripletsApproxWindow(const std::vector<HitVec>&, const BinInfoLayerMap&, const PairIdxVec&, TripletIdxVec&);
void intersectThirdLayer(const float, const float, const float, const float, float&, float&);
void filterHitTripletsBySecondLayerZResidual(const std::vector<HitVec>&, const TripletIdxVec&, TripletIdxVec&);
void filterHitTripletsByRZChi2(const std::vector<HitVec>&, const TripletIdxVec&, TripletIdxVec&);
void buildSeedsFromTriplets(const std::vector<HitVec>&, const TripletIdxVec&, TrackVec&, TrackExtraVec&, Event&);

#endif
