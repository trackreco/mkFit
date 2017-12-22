#ifndef _event_
#define _event_

#include "Track.h"
#include "Validation.h"
#include "Geometry.h"
#include "BinInfoUtils.h"
#include "Config.h"

#include <mutex>

struct DataFile;

class Event
{
public:
  explicit Event(int evtID);
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

  int  nextMCHitID() { return mcHitIDCounter_++; }

  void write_out(DataFile &data_file);
  void read_in  (DataFile &data_file, FILE *in_fp=0);

  void kludge_cms_hit_errors();

  int  use_seeds_from_cmsswtracks(); //special mode --> use only seeds which generated cmssw reco track
  int  clean_cms_simtracks();
  int  clean_cms_seedtracks(); //operates on seedTracks_; returns the number of cleaned seeds
  int  clean_cms_seedtracks_badlabel(); //operates on seedTracks_, removes those with label == -1;
  void relabel_bad_seedtracks();
  void print_tracks(const TrackVec& tracks, bool print_hits) const;

  const Geometry& geom_;
  Validation& validation_;

private:
  int  evtID_;

public:
  int threads_;
  std::mutex       mcGatherMutex_;
  std::atomic<int> mcHitIDCounter_;
  std::vector<HitVec> layerHits_;
  MCHitInfoVec simHitsInfo_;

  TrackVec simTracks_, seedTracks_, candidateTracks_, fitTracks_;
  TrackVec cmsswTracks_;
  // validation sets these, so needs to be mutable
  mutable TrackExtraVec simTracksExtra_, seedTracksExtra_, candidateTracksExtra_, fitTracksExtra_;
  mutable TrackExtraVec cmsswTracksExtra_;

  // XXXXMT: Preliminary. Separators into seed/candidate arrays.
  // This will have to be extended for multi-pass tracking.
  int seedEtaSeparators_[5];
  int seedMinLastLayer_[5];
  int seedMaxLastLayer_[5];

  TSVec simTrackStates_;
  static std::mutex printmutex;
};

typedef std::vector<Event> EventVec;


struct DataFileHeader
{
  int f_magic          = 0xBEEF;
  int f_format_version = 3;
  int f_sizeof_track   = sizeof(Track);
  int f_n_max_trk_hits = Config::nMaxTrkHits;
  int f_n_layers       = -1;
  int f_n_events       = -1;

  int f_extra_sections = 0;

  DataFileHeader()
  {
    f_n_layers = Config::nTotalLayers;
  }
};

struct DataFile
{
  enum ExtraSection
  {
    ES_SimTrackStates = 0x1,
    ES_Seeds          = 0x2,
    ES_CmsswTracks    = 0x4
  };

  FILE *f_fp  =  0;
  long  f_pos =  sizeof(DataFileHeader);

  DataFileHeader f_header;

  std::mutex     f_next_ev_mutex;

  // ----------------------------------------------------------------

  bool HasSimTrackStates() const { return f_header.f_extra_sections & ES_SimTrackStates; }
  bool HasSeeds()          const { return f_header.f_extra_sections & ES_Seeds; }
  bool HasCmsswTracks()    const { return f_header.f_extra_sections & ES_CmsswTracks; }

  int  OpenRead (const std::string& fname, bool set_n_layers = false);
  void OpenWrite(const std::string& fname, int nev, int extra_sections=0);

  int  AdvancePosToNextEvent(FILE *fp);

  void SkipNEvents(int n_to_skip);

  void Close();
  void CloseWrite(int n_written); //override nevents in the header and close
};

#endif
