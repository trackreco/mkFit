#ifndef MkFitter_h
#define MkFitter_h

#include "Matrix.h"
#include "KalmanUtils.h"


class MkFitter
{
  MPlexLS Err[2];
  MPlexLV Par[2];

  MPlexQI Chg;

  MPlexHS msErr[MAX_HITS];
  MPlexHV msPar[MAX_HITS];

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
  void InputTracksOnly   (std::vector<Track>& tracks, int beg, int end);
  void FitTracks();
  void OutputFittedTracks(std::vector<Track>& tracks, int beg, int end);

};

#endif
