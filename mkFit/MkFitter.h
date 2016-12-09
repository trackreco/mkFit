#ifndef MkFitter_h
#define MkFitter_h

#include "Event.h"
#include "Matrix.h"
#include "KalmanUtils.h"

#include "HitStructures.h"
#include "BinInfoUtils.h"

#if USE_CUDA
#include "FitterCU.h"
#include "HitStructuresCU.h"
#endif

//#define DEBUG 1

//#define USE_BOHS

class CandCloner;

const int MPlexHitIdxMax = 16;
typedef Matriplex::Matriplex<int, MPlexHitIdxMax, 1, NN> MPlexHitIdx;

struct MkFitter
{
  MPlexLS Err[2];
  MPlexLV Par[2];

  MPlexQI Chg;

  MPlexQF Chi2;

  MPlexHS msErr[Config::nLayers];
  MPlexHV msPar[Config::nLayers];

  MPlexQI Label;  //this is the seed index in global seed vector (for MC truth match)
  MPlexQI SeedIdx;//this is the seed index in local thread (for bookkeeping at thread level)
  MPlexQI CandIdx;//this is the candidate index for the given seed (for bookkeeping of clone engine)
  MPlexQI HitsIdx[Config::nLayers];

  // Hold hit indices to explore at current layer.
  MPlexQI     XHitSize;
  MPlexHitIdx XHitArr;

  // Indices into Err and Par arrays.
  // Thought I'll have to flip between them ...
  const int iC = 0; // current
  const int iP = 1; // propagated

  int Nhits;

public:
  MkFitter() : Nhits(0)
  {}
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

  float getPar(int itrack, int i, int par) const { return Par[i].ConstAt(itrack, 0, par); }
  void  SetNhits(int newnhits) { Nhits=newnhits; }

  int countValidHits  (int itrack, int end_hit) const;
  int countInvalidHits(int itrack, int end_hit) const;
  int countValidHits  (int itrack) const { return countValidHits  (itrack, Nhits); }
  int countInvalidHits(int itrack) const { return countInvalidHits(itrack, Nhits); }

  void InputTracksAndHits(const std::vector<Track>& tracks, const std::vector<HitVec>& layerHits, int beg, int end);
  void InputTracksAndHits(const std::vector<Track>& tracks, const std::vector<LayerOfHits>& layerHits, int beg, int end);
  void SlurpInTracksAndHits(const std::vector<Track>&  tracks, const std::vector<HitVec>& layerHits, int beg, int end);
  void InputTracksAndHitIdx(const std::vector<Track>& tracks,
                            int beg, int end, bool inputProp);
  void InputTracksAndHitIdx(const std::vector<std::vector<Track> >& tracks, const std::vector<std::pair<int,int> >& idxs,
                            int beg, int end, bool inputProp);
  void InputSeedsTracksAndHits(const std::vector<Track>& seeds, const std::vector<Track>& tracks, const std::vector<HitVec>& layerHits, int beg, int end);
  void ConformalFitTracks(bool fitting, int beg, int end);
  void FitTracks(const int N_proc, const Event * ev, const bool useParamBfield = false);
  void FitTracksTestEndcap(const int N_proc, const Event* ev, const bool useParamBfield = false);

  void CollectFitValidation(const int hi, const int N_proc, const Event * ev) const;

  void OutputTracks(std::vector<Track>& tracks, int beg, int end, int iCP) const;

  void OutputFittedTracks(std::vector<Track>& tracks, int beg, int end) const
  { return OutputTracks(tracks,beg,end,iC); }

  void OutputPropagatedTracks(std::vector<Track>& tracks, int beg, int end) const
  { return OutputTracks(tracks,beg,end,iP); }

  void OutputFittedTracksAndHitIdx(std::vector<Track>& tracks, int beg, int end, bool outputProp) const;

  void PropagateTracksToR(float R, const int N_proc);

  void PropagateTracksToZ(float Z, const int N_proc);

  void SelectHitIndices(const LayerOfHits &layer_of_hits, const int N_proc, bool dump=false);
  void SelectHitIndicesEndcap(const LayerOfHits &layer_of_hits, const int N_proc, bool dump=false);

  void AddBestHit      (const LayerOfHits &layer_of_hits, const int N_proc);
  void AddBestHitEndcap(const LayerOfHits &layer_of_hits, const int N_proc);

  void FindCandidates(const LayerOfHits &layer_of_hits, std::vector<std::vector<Track> >& tmp_candidates,
		      const int offset, const int N_proc);
  void FindCandidatesEndcap(const LayerOfHits &layer_of_hits, std::vector<std::vector<Track> >& tmp_candidates,
			    const int offset, const int N_proc);
  // ================================================================
  // Methods used with clone engine
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
  void FindCandidatesMinimizeCopy(const LayerOfHits &layer_of_hits, CandCloner& cloner,
                                  const int offset, const int N_proc);
  void FindCandidatesMinimizeCopyEndcap(const LayerOfHits &layer_of_hits, CandCloner& cloner,
                                        const int offset, const int N_proc);

  //version of input tracks using IdxChi2List
  void InputTracksAndHitIdx(const std::vector<std::vector<Track> >& tracks,
                            const std::vector<std::pair<int,IdxChi2List> >& idxs,
                            int beg, int end, bool inputProp = false);

  void UpdateWithLastHit(const LayerOfHits &layer_of_hits, int N_proc);
  void UpdateWithLastHitEndcap(const LayerOfHits &layer_of_hits, int N_proc);

  void CopyOutParErr(std::vector<std::vector<Track> >& seed_cand_vec,
                     int N_proc, bool outputProp) const;
};

#endif
