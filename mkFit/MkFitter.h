#ifndef MkFitter_h
#define MkFitter_h

#include "Matrix.h"
#include "KalmanUtils.h"

//#define DEBUG 1

class MkFitter
{
  MPlexLS Err[2];
  MPlexLV Par[2];

  MPlexQI Chg;

  MPlexQF Chi2;

  MPlexHS msErr[MAX_HITS];
  MPlexHV msPar[MAX_HITS];

  MPlexQI SeedIdx;
  MPlexQI HitsIdx[MAX_HITS];

  // Indices into Err and Par arrays.
  // Thought I'll have to flip between them ...
  const int iC = 0; // current
  const int iP = 1; // propagated

  int Nhits;

public:
  MkFitter(int n_hits) : Nhits(n_hits)
  {
    // XXXX Eventually dynamically allocate measurement arrays.
    // XXXX std::vector is no good, not aligned!
    // XXXX Hmmh, should really copy them in layer by layer.
  }

  // Copy-in timing tests.
  MPlexLS& GetErr0() { return Err[0]; }
  MPlexLV& GetPar0() { return Par[0]; }

  void CheckAlignment();

  void PrintPt(int idx);

  void InputTracksAndHits(std::vector<Track>& tracks, int beg, int end);
  void InputTracksAndHitIdx(std::vector<std::vector<Track> >& tracks, std::vector<std::pair<int,int> >& idxs, int beg, int end);
  void InputTracksOnly   (std::vector<Track>& tracks, int beg, int end);
  void InputHitsOnly(std::vector<Hit>& hits, int beg, int end);
  void FitTracks();

  void OutputTracks(std::vector<Track>& tracks, int beg, int end, int iCP);
  void OutputFittedTracks(std::vector<Track>& tracks, int beg, int end) {
    return OutputTracks(tracks,beg,end,iC);
  }
  void OutputPropagatedTracks(std::vector<Track>& tracks, int beg, int end){
    return OutputTracks(tracks,beg,end,iP);
  }

  void OutputFittedTracksAndHits(std::vector<Track>& tracks, int beg, int end);

  void PropagateTracksToR(float R);

  void AddBestHit(std::vector<Hit>& lay_hits, int beg, int end);

  typedef std::pair<unsigned int,unsigned int> BinInfo;
  void GetHitRange(std::vector<BinInfo>& segmentMapLay_, int beg, int end, const float& etaDet, int& firstHit, int& lastHit);

  void FindCandidates(std::vector<Hit>& lay_hits, int firstHit, int lastHit, int beg, int end, std::vector<std::pair<Track,int> >& reccands_tmp);

  void SetNhits(int newnhits) { Nhits=newnhits; }

  int countInvalidHits(int itrack);

};

#endif
