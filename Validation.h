#ifndef _validation_
#define _validation_
#include "Track.h"
//#include <map>

class Validation {
public:
  virtual void fillBuildTree(unsigned int, unsigned int, unsigned int) {}
  virtual void fillEffTree(TrackVec&, TrackVec&, TrackVec&, TrackVec&) {}
  virtual void saveTTrees() {}
};

#endif
