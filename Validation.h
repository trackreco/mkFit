#ifndef _validation_
#define _validation_
#include "Track.h"

class Validation {
public:
  virtual void fillSimHists(TrackVec& evt_seeds) {}
  virtual void fillCandidateHists(TrackVec& evt_track_candidates) {}
  virtual void fillBuildHists(unsigned int, unsigned int, unsigned int) {}
  virtual void fillFitStateHists(TrackState&, TrackState&) {}
  virtual void fillFitHitHists(MeasurementState&, MeasurementState&, TrackState&, TrackState&) {}
  virtual void fillFitTrackHists(TrackState&, TrackState&) {}
  virtual void saveHists() {}
  virtual void deleteHists() {}
};

#endif
