#ifndef _validation_
#define _validation_
#include "Track.h"

class Validation {
public:
  virtual void fillSimHists(TrackVec& evt_sim_tracks) {}
  virtual void fillSeedHists(std::vector<HitVec> &, std::vector<HitVec> &, TrackVec&, HitVec &, HitVec &, std::vector<float> &, std::vector<float> &, TrackVec &, TrackVec &, TrackVec &) {}
  virtual void fillCandidateHists(TrackVec& evt_cand_tracks) {}
  virtual void fillBuildHists(unsigned int, unsigned int, unsigned int) {}
  virtual void fillAssociationHists(TrackVec& , TrackVec& ) {}
  virtual void fillFitStateHists(TrackState&, TrackState&) {}
  virtual void fillFitHitHists(unsigned int, HitVec&, MeasurementState&, TrackState&, TrackState&) {}
  virtual void fillFitTrackHists(TrackState&, TrackState&) {}
  virtual void saveHists() {}
  virtual void deleteHists() {}
};

#endif
