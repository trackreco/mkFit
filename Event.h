#ifndef _event_
#define _event_

#include "Track.h"
#include "Validation.h"
#include "Geometry.h"
#include "BinInfoUtils.h"
#include "Config.h"

struct HitID {
  HitID() : layer(-1), index(-1) {}
  HitID(short l, short i) : layer(l), index(i) {}
  short layer;
  short index;
};
typedef std::vector<HitID> HitIDVec;

class Event {
public:
  Event(const Geometry& g, Validation& v, unsigned int evtID, int threads = 1);
  void Simulate();
  void Segment();
  void Seed();
  void Find();
  void Fit();
  void Validate(const unsigned int);
  
  unsigned int evtID() const {return evtID_;}
  void resetLayerHitMap(bool resetSimHits);

  const Geometry& geom_;
  Validation& validation_;
 private:
  unsigned int evtID_;
 public:
  int threads_;
  std::vector<HitVec> layerHits_;
  MCHitInfoVec simHitsInfo_;
  HitIDVec layerHitMap_; // indexed same as simHitsInfo_, maps to layer & hit

  TrackVec simTracks_, seedTracks_, candidateTracks_, fitTracks_;

  // phi-eta partitioning map: vector of vector of vectors of std::pairs. 
  // vec[nLayers][nEtaBins][nPhiBins]
  BinInfoMap segmentMap_;
};

typedef std::vector<Event> EventVec;

#endif
