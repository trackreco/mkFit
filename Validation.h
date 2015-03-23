#ifndef _validation_
#define _validation_
#include "Track.h"

class Validation {
public:
  virtual void fillSimHists(const TrackVec&) {}
  virtual void fillSeedHists(std::vector<HitVec> &, std::vector<HitVec> &, TrackVec&, HitVec &, HitVec &, std::vector<float> &, std::vector<float> &, TrackVec &, TrackVec &, TrackVec &) {}
  virtual void fillCandidateHists(const TrackVec&) {}
  virtual void fillBuildHists(unsigned int, unsigned int, unsigned int) {}
  virtual void fillAssociationHists(const TrackVec& evt_track_candidates, const TrackVec& evt_sim_tracks) {}
  virtual void fillFitStateHists(const TrackState&, const TrackState&) {}
  virtual void fillFitHitHists(unsigned int, const HitVec&, const MeasurementState&, const TrackState&, const TrackState&) {}
  virtual void fillFitTrackHists(const TrackState&, const TrackState&) {}
  virtual void saveHists() {}
  virtual void deleteHists() {}
};

#endif
