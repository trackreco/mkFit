#ifndef _event_
#define _event_

#include "Track.h"
#include "Validation.h"
#include "Geometry.h"
#include "BinInfoUtils.h"
#include "Config.h"

#include <mutex>

class Event
{
public:
  Event(const Geometry& g, Validation& v, int evtID, int threads = 1);
  void Reset(int evtID);
  void RemapHits(TrackVec & tracks);
  void Simulate();
  void Segment(BinInfoMap & segmentMap);
  void Seed(const BinInfoMap & segmentMap);
  void Find(const BinInfoMap & segmentMap);
  void Fit();
  void Validate();
  void PrintStats(const TrackVec&, TrackExtraVec&);
  
  int  evtID() const {return evtID_;}
  void resetLayerHitMap(bool resetSimHits);

  int nextMCHitID() { return mcHitIDCounter_++; }

  void write_out(FILE *fp);
  void read_in(FILE *fp, int version = Config::FileVersion);

  void clean_cms_simtracks();
  void print_tracks(const TrackVec& tracks, bool print_hits) const;

  const Geometry& geom_;
  Validation& validation_;

private:
  int evtID_;

public:
  int threads_;
  std::mutex       mcGatherMutex_;
  std::atomic<int> mcHitIDCounter_;
  std::vector<HitVec> layerHits_;
  MCHitInfoVec simHitsInfo_;

  TrackVec simTracks_, seedTracks_, candidateTracks_, fitTracks_;
  // validation sets these, so needs to be mutable
  mutable TrackExtraVec simTracksExtra_, seedTracksExtra_, candidateTracksExtra_, fitTracksExtra_;

  // XXXXMT: Preliminary. Separators into seed/candidate arrays.
  // This will have to be extended for multi-pass tracking.
  int seedEtaSeparators_[5];
  int seedMinLastLayer_[5];

  TSVec simTrackStates_;
  static std::mutex printmutex;
};

typedef std::vector<Event> EventVec;

#endif
