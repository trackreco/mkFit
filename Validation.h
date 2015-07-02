#ifndef _validation_
#define _validation_
#include "Track.h"

class Validation {
public:
  virtual void resetValidationMaps() {}

  virtual void collectSimTkTSVecMapInfo(const unsigned int, const TSVec&) {}

  virtual void collectSeedTkCFMapInfo(const unsigned int, const TrackState&) {}
  virtual void collectSeedTkTSLayerPairVecMapInfo(const unsigned int, const TSLayerPairVec&) {}

  virtual void collectFitTkCFMapInfo(const unsigned int, const TrackState&) {}
  virtual void collectFitTkTSLayerPairVecMapInfo(const unsigned int, const TSLayerPairVec&) {}

  virtual void fillBuildTree(const unsigned int, const unsigned int, const std::vector<unsigned int>&, const std::vector<unsigned int>&, const std::vector<unsigned int>&) {}

  virtual void makeSimTkToRecoTksMaps(TrackVec&, TrackVec&, TrackVec&) {}
  virtual void fillEffTree(const TrackVec&, const unsigned int) {}

  virtual void makeSeedTkToRecoTkMaps(const TrackVec&, const TrackVec&) {}
  virtual void fillFakeRateTree(const TrackVec&, const unsigned int) {}

  virtual void fillConfigTree(const std::vector<double> &) {}

  virtual void saveTTrees() {}
};

#endif
