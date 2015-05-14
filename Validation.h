#ifndef _validation_
#define _validation_
#include "Track.h"

class Validation {
public:
  virtual void fillBuildTree(const unsigned int, const unsigned int, const unsigned int) {}
  virtual void makeSimToTkMaps(TrackVec&, TrackVec&, TrackVec&) {}
  virtual void fillEffTree(const TrackVec&, const unsigned int &) {}
  virtual void fillFakeTrees(const unsigned int &) {}
  virtual void fillConfigTree(const unsigned int &, const unsigned int &, const std::vector<double> &) {}
  virtual void saveTTrees() {}
};

#endif
