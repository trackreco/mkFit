#ifndef MkFitter_h
#define MkFitter_h

#include "Event.h"
#include "Matrix.h"
#include "KalmanUtils.h"

#include "HitStructures.h"
#include "BinInfoUtils.h"

//#define DEBUG 1

class CandCloner;

class MkFitter
{
  MPlexLS Err[2];
  MPlexLV Par[2];

  MPlexQI Chg;

  MPlexQF Chi2;

  MPlexHS msErr[MAX_HITS];
  MPlexHV msPar[MAX_HITS];

  MPlexQI Label;  //this is the seed index in global seed vector (for MC truth match)
  MPlexQI SeedIdx;//this is the seed index in local thread (for bookkeeping at thread level)
  MPlexQI CandIdx;//this is the candidate index for the given seed (for bookkeeping of clone engine)
  MPlexQI HitsIdx[MAX_HITS];

  // Hold hit indices to explore at current layer.
  MPlexQI XHitPos;   // Should be pos (so can move it forward)
  MPlexQI XHitSize;  // Should be size (so i can reduce it)
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

  //int getXHitBegin(int arg0,int arg1,int arg2) { return XHitBegin.At(arg0, arg1, arg2); }
  //int getXHitEnd  (int arg0,int arg1,int arg2) { return XHitEnd  .At(arg0, arg1, arg2); }

  void InputTracksAndHits(std::vector<Track>& tracks, std::vector<HitVec>& layerHits, int beg, int end);
  void InputTracksAndHitIdx(std::vector<Track>& tracks,
                            int beg, int end, bool inputProp);
  void InputTracksAndHitIdx(std::vector<std::vector<Track> >& tracks, std::vector<std::pair<int,int> >& idxs,
                            int beg, int end, bool inputProp);
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

  void OutputFittedTracksAndHitIdx(std::vector<Track>& tracks, int beg, int end, bool outputProp);

  void PropagateTracksToR(float R, const int N_proc);

  void AddBestHit(std::vector<Hit>& lay_hits, int firstHit, int lastHit, int beg, int end);

  void GetHitRange(std::vector<std::vector<BinInfo> >& segmentMapLay_, int beg, int end,
                   int& firstHit, int& lastHit);

  void FindCandidates(std::vector<Hit>& lay_hits, int firstHit, int lastHit, int beg, int end, std::vector<std::vector<Track> >& tmp_candidates, int offset);

  void SetNhits(int newnhits) { Nhits=newnhits; }

  int countValidHits  (int itrack, int end_hit);
  int countInvalidHits(int itrack, int end_hit);
  int countValidHits  (int itrack) { return countValidHits  (itrack, Nhits); }
  int countInvalidHits(int itrack) { return countInvalidHits(itrack, Nhits); }

  float getPar(int itrack, int i, int par) { return Par[i].ConstAt(itrack, 0, par); }


  // ================================================================
  // MT methods
  // ================================================================

  void SelectHitRanges(BunchOfHits &bunch_of_hits, const int N_proc);
  void AddBestHit     (BunchOfHits &bunch_of_hits);

  void FindCandidates(BunchOfHits &bunch_of_hits, std::vector<std::vector<Track> >& tmp_candidates,
                      const int offset, const int N_proc);

  // ================================================================
  // Methods to be used with clone engine
  // ================================================================
  //minimal set of information for bookkeping
  struct IdxChi2List
  {
    int   trkIdx;//candidate index
    int   hitIdx;//hit index
    int   nhits; //number of hits (used for sorting)
    float chi2;//total chi2 (used for sorting)
  };
  //version of find candidates that does not cloning, just fills the IdxChi2List as output (to be then read by the clone engine)
  void FindCandidatesMinimizeCopy(BunchOfHits &bunch_of_hits, CandCloner& cloner,
                                  const int offset, const int N_proc);

  //version of input tracks using IdxChi2List
  void InputTracksAndHitIdx(std::vector<std::vector<Track> >& tracks,
                            std::vector<std::pair<int,IdxChi2List> >& idxs,
                            int beg, int end, bool inputProp = false);

  //method used by the clone engine to do the actual cloning on the predefined candidate+hit
  void UpdateWithHit(BunchOfHits &bunch_of_hits,
		     std::vector<std::pair<int,IdxChi2List> >& idxs,
		     std::vector<std::vector<Track> >& cands_for_next_lay,
		     int offset, int beg, int end);

  //method used by the clone engine to do the actual cloning on the predefined candidate+hit
  void UpdateWithHit(BunchOfHits &bunch_of_hits,
                     std::vector<std::pair<int,IdxChi2List> >& idxs,
                     int beg, int end);

  //method used by the clone engine to do the actual cloning on the predefined candidate+hit
  void CopyOutClone(std::vector<std::pair<int,IdxChi2List> >& idxs,
		    std::vector<std::vector<Track> >& cands_for_next_lay,
		    int offset, int beg, int end, bool outputProp = false);
};

#endif
