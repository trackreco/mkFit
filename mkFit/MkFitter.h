#ifndef MkFitter_h
#define MkFitter_h

#include "Event.h"
#include "Matrix.h"
#include "KalmanUtils.h"

#include "HitStructures.h"

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

  // Hold hit indices to explore at current layer.
  MPlexQI XHitBegin; // Should be pos (so can move it forward)
  MPlexQI XHitEnd;   // Should be size (so i can reduce it)
                     // with pos size could even handle phi wrap! XXXX do it
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

  int getXHitBegin(int arg0,int arg1,int arg2) { return XHitBegin.At(arg0, arg1, arg2); }
  int getXHitEnd  (int arg0,int arg1,int arg2) { return XHitEnd  .At(arg0, arg1, arg2); }

  void InputTracksAndHits(std::vector<Track>& tracks, int beg, int end);
  void InputTracksAndHitIdx(std::vector<Track>& tracks, int beg, int end);
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
  void OutputFittedTracksAndHitIdx(std::vector<Track>& tracks, int beg, int end);

  void PropagateTracksToR(float R);

  void AddBestHit(std::vector<Hit>& lay_hits, int firstHit, int lastHit, int beg, int end);

  void GetHitRange(std::vector<std::vector<BinInfo> >& segmentMapLay_, int beg, int end,
                   const float etaDet, int& firstHit, int& lastHit);

  void FindCandidates(std::vector<Hit>& lay_hits, int firstHit, int lastHit, int beg, int end, std::vector<std::vector<Track> >& tmp_candidates, int offset);

  void SetNhits(int newnhits) { Nhits=newnhits; }

  int countInvalidHits(int itrack);

  float getPar(int itrack, int i, int par) { return Par[i].ConstAt(itrack, 0, par); }

  inline float normalizedPhi(float phi) {
    return std::fmod(phi, (float) M_PI);
  }


  // ================================================================
  // MT methods
  // ================================================================

  void SelectHitRanges(BunchOfHits &bunch_of_hits);
  void AddBestHit     (BunchOfHits &bunch_of_hits);
};

#endif
