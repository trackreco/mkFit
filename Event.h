#ifndef _event_
#define _event_

#include "Track.h"
#include "Validation.h"
#include "Geometry.h"
#include "BinInfoUtils.h"

namespace Config {
  static constexpr const unsigned int nlayers_per_seed = 3;
  static constexpr const unsigned int maxCand = 10;
  static constexpr const float chi2Cut = 15.;
  static constexpr const float nSigma = 3.;
  static constexpr const float minDPhi = 0.;
};

class Event {
public:
  Event(const Geometry& g, Validation& v, int threads = 1);
  void Simulate(unsigned int nTracks);
  void Segment();
  void Seed();
  void Find();
  void Fit();
  void ValidateHighLevel(const unsigned int &);

  const Geometry& geom_;
  Validation& validation_;
  std::vector<HitVec> layerHits_;
  TrackVec simTracks_, seedTracks_, candidateTracks_, fitTracks_;
  int threads_;

  BinInfoMap segmentMap_;
};

typedef std::vector<Event> EventVec;

#endif
