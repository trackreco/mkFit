#ifndef _event_
#define _event_

#include "Track.h"
#include "Validation.h"
#include "Geometry.h"

typedef std::pair<unsigned int,unsigned int> BinInfo;

/*
struct LayEtaPhiInfo{
  BinInfo lay_eta_hit_idx;
  std::vector<BinInfo> vec_lay_phi_hit_idx;
};
*/

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
  //  TrackVec simTracks_, seedTracks_, candidateTracks_, associatedTracks_RD_, associatedTracks_SD_, fitTracks_;

  //these matrices are dummy and can be optimized without multiplying by zero all the world...
  SMatrix36 projMatrix36_;
  SMatrix63 projMatrix36T_;

  // phi partitioning map: Vector of vectors of std::pairs. 
  // A vector of maps, although vector is fixed to layer, so really array of maps,
  // where maps are phi bins and the number of hits in those phi bins.
  // First is first hit index in bin, second is size of this bin
  //  std::vector<std::vector<BinInfo> > lay_phi_hit_idx_;
  //  std::vector<std::vector<LayEtaPhiInfo> > lay_eta_phi_hit_idx_;

  std::vector<std::vector<std::vector<BinInfo> > > lay_eta_phi_hit_idx_;
};

typedef std::vector<Event> EventVec;

#endif
