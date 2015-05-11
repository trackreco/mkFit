#ifndef _validation_
#define _validation_
#include "Track.h"

class Validation {
public:
  virtual void fillBuildTree(const unsigned int, const unsigned int, const unsigned int) {}
  virtual void makeSimToTkMaps(TrackVec&, TrackVec&, TrackVec&) {}
  virtual void fillEffTree(const TrackVec&, const TrackVec&, const TrackVec&, const TrackVec&, const unsigned int &) {}
  virtual void saveTTrees() {}
};

#endif
