#ifndef _event_
#define _event_

#include "Track.h"
#include "Validation.h"
#include "Geometry.h"

typedef std::pair<unsigned int,unsigned int> BinInfo;

class Event {
public:
  Event(Geometry& g, Validation& v);
  void Simulate(unsigned int nTracks);
  void Segment();
  void Seed();
  void Find();
  void Fit();

  Geometry& geom_;
  Validation& validation_;
  std::vector<HitVec> layerHits_;
  TrackVec simTracks_, seedTracks_, candidateTracks_;

  //these matrices are dummy and can be optimized without multiplying by zero all the world...
  SMatrix36 projMatrix36_;
  SMatrix63 projMatrix36T_;

  // phi-eta partitioning map: vector of vector of vectors of std::pairs. 
  // vec[nLayers][nEtaBins][nPhiBins]
  std::vector<std::vector<std::vector<BinInfo> > > lay_eta_phi_hit_idx_;

};

typedef std::vector<Event> EventVec;

inline float normalizedPhi(float phi) {
  static float const TWO_PI = M_PI * 2;
  while ( phi < -M_PI ) phi += TWO_PI;
  while ( phi >  M_PI ) phi -= TWO_PI;
  return phi;
}

inline float normalizedEta(float eta) {
  static float const ETA_DET = 2.0;

  if (eta < -ETA_DET ) eta = -ETA_DET+.00001;
  if (eta >  ETA_DET ) eta =  ETA_DET-.00001;
  return eta;
}

#endif
