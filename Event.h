#ifndef _event_
#define _event_

#include "Track.h"
#include "Validation.h"
#include "Geometry.h"
#include "BinInfoUtils.h"
#include "Config.h"

class Event {
public:
  Event(const Geometry& g, Validation& v, unsigned int evtID, int threads = 1);
  void Simulate();
  void Segment();
  void Seed();
  void Find();
  void Fit();
  void Validate(const unsigned int);
  
  const unsigned int evtID() const {return evtID_;}
  
  const Geometry& geom_;
  Validation& validation_;
 private:
  unsigned int evtID_;
 public:
  int threads_;
  std::vector<HitVec> layerHits_;
  TrackVec simTracks_, seedTracks_, candidateTracks_, fitTracks_;

  // phi-eta partitioning map: vector of vector of vectors of std::pairs. 
  // vec[nLayers][nEtaBins][nPhiBins]
  BinInfoMap segmentMap_;
};

typedef std::vector<Event> EventVec;

#endif
