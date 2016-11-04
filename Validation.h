#ifndef _validation_
#define _validation_
#include "Track.h"
#include "BinInfoUtils.h"

class Event;

// Fit Validation objects -- mplex only
struct FitVal
{
public:
  FitVal() {}
  FitVal(float ppz, float eppz, float ppphi, float eppphi, 
	 float upt, float eupt, float umphi, float eumphi, float umeta, float eumeta) :
  ppz(ppz), eppz(eppz), ppphi(ppphi), eppphi(eppphi), upt(upt), eupt(eupt), umphi(umphi), eumphi(eumphi), umeta(umeta), eumeta(eumeta) {}

  // first p or u = propagated or updated
  // middle: p or m/nothing = position or momentum
  // begining: e = error (already sqrt)
  float ppz, eppz, ppphi, eppphi;
  float upt, eupt, umphi, eumphi, umeta, eumeta;
};

class Validation {
public:
  virtual void alignTrackExtra(TrackVec&, TrackExtraVec&) {}

  virtual void resetValidationMaps() {}
  virtual void resetDebugVectors() {}

  virtual void collectSimTkTSVecMapInfo(int, const TSVec&) {}
  virtual void collectSeedTkCFMapInfo(int, const TrackState&) {}
  virtual void collectSeedTkTSLayerPairVecMapInfo(int, const TSLayerPairVec&) {}
  virtual void collectBranchingInfo(int, int, float, float, int, float, int,
				    int, const std::vector<int>&, const std::vector<int>&) {}
  virtual void collectFitTkCFMapInfo(int, const TrackState&) {}
  virtual void collectFitTkTSLayerPairVecMapInfo(int, const TSLayerPairVec&) {}

  virtual void collectFitInfo(const FitVal&, int, int) {}

  virtual void collectPropTSLayerVecInfo(int, const TrackState&) {} // exclusively for debugtree
  virtual void collectChi2LayerVecInfo(int, float) {} // exclusively for debugtree
  virtual void collectUpTSLayerVecInfo(int, const TrackState&) {} // exclusively for debugtree

  virtual void setTrackExtras(Event& ev) {}
  virtual void makeSimTkToRecoTksMaps(Event&) {}
  virtual void makeSeedTkToRecoTkMaps(Event&) {}

  virtual void fillSeedInfoTree(const TripletIdxVec&, const Event&) {}
  virtual void fillSeedTree(const TripletIdxVec&, const TripletIdxVec&, const Event&) {}
  virtual void fillDebugTree(const Event&) {}
  virtual void fillSegmentTree(const BinInfoMap&, int) {}
  virtual void fillBranchTree(int) {}
  virtual void fillEfficiencyTree(const Event&) {}
  virtual void fillFakeRateTree(const Event&) {}
  virtual void fillGeometryTree(const Event&) {}
  virtual void fillConformalTree(const Event&) {}
  virtual void fillConfigTree() {}
  virtual void fillTimeTree(const std::vector<double> &) {}
  virtual void fillFitTree(const Event&) {}

  virtual void saveTTrees() {}
};

#endif
