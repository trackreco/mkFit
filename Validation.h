#ifndef _validation_
#define _validation_
#include "Track.h"
#include "BinInfoUtils.h"

class Event;

class Validation {
public:
  virtual void resetValidationMaps() {}

  virtual void collectSimTkTSVecMapInfo(unsigned int, const TSVec&) {}

  virtual void collectSeedTkCFMapInfo(unsigned int, const TrackState&) {}
  virtual void collectSeedTkTSLayerPairVecMapInfo(unsigned int, const TSLayerPairVec&) {}

  virtual void collectBranchingInfo(unsigned int, unsigned int, float, float, unsigned int, float, unsigned int, unsigned int, const std::vector<unsigned int>&, const std::vector<unsigned int>&) {}

  virtual void collectFitTkCFMapInfo(unsigned int, const TrackState&) {}
  virtual void collectFitTkTSLayerPairVecMapInfo(unsigned int, const TSLayerPairVec&) {}

  virtual void fillSegmentTree(const BinInfoMap&, unsigned int) {}

  virtual void fillBranchTree(unsigned int) {}

  virtual void makeSimTkToRecoTksMaps(const Event&) {}
  virtual void fillEffTree(const TrackVec&, unsigned int) {}

  virtual void makeSeedTkToRecoTkMaps(const TrackVec&, const TrackVec&) {}
  virtual void fillFakeRateTree(const TrackVec&, unsigned int) {}

  virtual void fillConfigTree(const std::vector<double> &) {}

  virtual void saveTTrees() {}
};

#endif
