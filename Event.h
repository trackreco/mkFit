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
  static constexpr const float maxDPhi = M_PI;
  static constexpr const float minDEta = 0.;
  static constexpr const float maxDEta = 1.0;
};

class Event {
public:
  Event(const Geometry& g, Validation& v, unsigned int evtID, int threads = 1);
  void Simulate(unsigned int nTracks);
  void Segment();
  void Seed();
  void Find();
  void Fit();
  void ValidateHighLevel(const unsigned int &);
  
  const unsigned int evtID() const {return evtID_;}
  
  const Geometry& geom_;
  Validation& validation_;
  std::vector<HitVec> layerHits_;
  TrackVec simTracks_, seedTracks_, candidateTracks_, fitTracks_;
  int threads_;

  // phi-eta partitioning map: vector of vector of vectors of std::pairs. 
  // vec[nLayers][nEtaBins][nPhiBins]
  BinInfoMap segmentMap_;
 private:
  unsigned int evtID_;
};

typedef std::vector<Event> EventVec;

#endif
