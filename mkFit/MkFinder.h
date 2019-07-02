#ifndef MkFinder_h
#define MkFinder_h

#include "MkBase.h"
#include "TrackerInfo.h"
#include "Track.h"

//#include "Event.h"

// Needed for TrackCand
#include "HitStructures.h"

namespace mkfit {

class CandCloner;
class CombCandidate;
class LayerOfHits;
class FindingFoos;

// For backward fit hack
class EventOfHits;
class EventOfCombCandidates;
class SteeringParams;

// NOTES from MkFitter ... where things were getting super messy.
//
// Would like to get rid of msErr, msPar arrays ... single one should do.
// Would also like to get rid of Nhits ... but this is a tougher one.
// Fitting is somewhat depending on this (not the new inter-slurp version).
// No need to have hit array as matriplex. Ideally, don't even want the full
// array if we could copy / clone it from the original candidate.
// Copy as needed or copy full size (compile time constant)? memcpy, std::copy 
// Do need per track size.
//
// Can also have nFoundHits ... then do not need countValid/InvalidHits().
//
// For now will prefix the new variables with 'tf' for track finding.
// Changes will go into all finding functions + import-with-hit-idcs.
//
// Actually ... am tempted to make MkFinder :)


// Define to get printouts about track and hit chi2.
// See also MkBuilder::BackwardFit() and MkBuilder::quality_store_tracks().

// #define DEBUG_BACKWARD_FIT


class MkFinder : public MkBase
{
public:

  static constexpr int MPlexHitIdxMax = 16;

  using MPlexHitIdx = Matriplex::Matriplex<int, MPlexHitIdxMax, 1, NN>;
  using MPlexQHoT   = Matriplex::Matriplex<HitOnTrack, 1, 1, NN>;

  //----------------------------------------------------------------------------

  MPlexQF    Chi2;
  MPlexQI    Label;   // seed index in global seed vector (for MC truth match)

  MPlexQI    NHits;
  MPlexQI    NFoundHits;
  HitOnTrack HoTArrs[NN][Config::nMaxTrkHits];

  MPlexQUI   SeedType; // seed range for ranking (0 = not set; 1 = high pT central seeds; 2 = low pT endcap seeds; 3 = all other seeds)
  MPlexQI    SeedIdx; // seed index in local thread (for bookkeeping at thread level)
  MPlexQI    CandIdx; // candidate index for the given seed (for bookkeeping of clone engine)

  // Additions / substitutions for TrackCand copy_in/out()
  MPlexQI    NMissingHits;   // sub: NHits, sort of
  MPlexQI    NInsideMinusOneHits;  // sub: before we copied all hit idcs and had a loop counting them
  MPlexQI    NTailMinusOneHits;  // sub: before we copied all hit idcs and had a loop counting them
  MPlexQI    LastHitCcIndex; // add: index of last hit in CombCand hit tree
  HitOnTrack LastHoT[NN];
  CombCandidate *CombCand[NN];
  // const TrackCand *TrkCand[NN]; // hmmh, could get all data through this guy ... but scattered

  // Hit indices into LayerOfHits to explore.
  WSR_Result  XWsrResult[NN]; // Could also merge it with XHitSize. Or use smaller arrays.
  MPlexQI     XHitSize;
  MPlexHitIdx XHitArr;

  // Hit errors / parameters for hit matching, update.
  MPlexHS    msErr;
  MPlexHV    msPar;

  // An idea: Do propagation to hit in FindTracksXYZZ functions.
  // Have some state / functions here that make this short to write.
  // This would simplify KalmanUtils (remove the propagate functions).
  // Track errors / parameters propagated to current hit.
  // MPlexLS    candErrAtCurrHit;
  // MPlexLV    candParAtCurrHit;

  //============================================================================

  MkFinder() {}

  //----------------------------------------------------------------------------

  void InputTracksAndHitIdx(const std::vector<Track>& tracks,
                            int beg, int end, bool inputProp);

  void InputTracksAndHitIdx(const std::vector<Track>& tracks,
                            const std::vector<int>  &   idxs,
                            int beg, int end, bool inputProp, int mp_offset);

  void InputTracksAndHitIdx(const std::vector<CombCandidate>& tracks,
                            const std::vector<std::pair<int,int>>& idxs,
                            int beg, int end, bool inputProp);

  void InputTracksAndHitIdx(const std::vector<CombCandidate>& tracks,
                            const std::vector<std::pair<int,IdxChi2List>>& idxs,
                            int beg, int end, bool inputProp);

  void OutputTracksAndHitIdx(std::vector<Track>& tracks,
                             int beg, int end, bool outputProp) const;

  void OutputTracksAndHitIdx(std::vector<Track>& tracks,
                             const std::vector<int>& idxs,
                             int beg, int end, bool outputProp) const;

  //----------------------------------------------------------------------------

  void getHitSelDynamicWindows(const LayerOfHits &layer_of_hits, const float track_pt, const float track_eta, float &min_dq, float &max_dphi);

  void SelectHitIndices(const LayerOfHits &layer_of_hits, const int N_proc);

  void AddBestHit(const LayerOfHits &layer_of_hits, const int N_proc,
                  const FindingFoos &fnd_foos);

  //----------------------------------------------------------------------------

  void FindCandidates(const LayerOfHits &layer_of_hits,
                      std::vector<std::vector<TrackCand>>& tmp_candidates,
		      const int offset, const int N_proc,
                      const FindingFoos &fnd_foos);

  //----------------------------------------------------------------------------

  void FindCandidatesCloneEngine(const LayerOfHits &layer_of_hits, CandCloner& cloner,
                                 const int offset, const int N_proc,
                                 const FindingFoos &fnd_foos);

  void UpdateWithLastHit(const LayerOfHits &layer_of_hits, int N_proc,
                         const FindingFoos &fnd_foos);

  void CopyOutParErr(std::vector<CombCandidate>& seed_cand_vec,
                     int N_proc, bool outputProp) const;

  //----------------------------------------------------------------------------
  // Backward fit hack

  int               CurHit[NN];
  const HitOnTrack *HoTArr[NN];

  void BkFitInputTracks (TrackVec& cands, int beg, int end);
  void BkFitOutputTracks(TrackVec& cands, int beg, int end);
  void BkFitInputTracks (EventOfCombCandidates& eocss, int beg, int end);
  void BkFitOutputTracks(EventOfCombCandidates& eocss, int beg, int end);

  void BkFitFitTracks(const EventOfHits& eventofhits, const SteeringParams& st_par,
                      const int N_proc, bool chiDebug = false);

  void BkFitPropTracksToPCA(const int N_proc);

  //----------------------------------------------------------------------------

private:

  void copy_in(const Track& trk, const int mslot, const int tslot)
  {
    Err[tslot].CopyIn(mslot, trk.errors().Array());
    Par[tslot].CopyIn(mslot, trk.parameters().Array());

    Chg  (mslot, 0, 0) = trk.charge();
    Chi2 (mslot, 0, 0) = trk.chi2();
    Label(mslot, 0, 0) = trk.label();

    NHits     (mslot, 0, 0) = trk.nTotalHits();
    NFoundHits(mslot, 0, 0) = trk.nFoundHits();
    std::copy(trk.BeginHitsOnTrack(), trk.EndHitsOnTrack(), HoTArrs[mslot]); 
  }

  void copy_out(Track& trk, const int mslot, const int tslot) const
  {
    Err[tslot].CopyOut(mslot, trk.errors_nc().Array());
    Par[tslot].CopyOut(mslot, trk.parameters_nc().Array());

    trk.setCharge(Chg  (mslot, 0, 0));
    trk.setChi2  (Chi2 (mslot, 0, 0));
    trk.setLabel (Label(mslot, 0, 0));

    trk.setNTotalHits(NHits     (mslot, 0, 0));
    trk.setNFoundHits(NFoundHits(mslot, 0, 0));
    std::copy(HoTArrs[mslot], & HoTArrs[mslot][NHits(mslot, 0, 0)], trk.BeginHitsOnTrack_nc());
  }

  void copy_in(const TrackCand& trk, const int mslot, const int tslot)
  {
    Err[tslot].CopyIn(mslot, trk.errors().Array());
    Par[tslot].CopyIn(mslot, trk.parameters().Array());

    Chg  (mslot, 0, 0) = trk.charge();
    Chi2 (mslot, 0, 0) = trk.chi2();
    Label(mslot, 0, 0) = trk.label();

    LastHitCcIndex(mslot, 0, 0) = trk.lastCcIndex();
    NFoundHits    (mslot, 0, 0) = trk.nFoundHits();
    NMissingHits  (mslot, 0, 0) = trk.nMissingHits();

    NInsideMinusOneHits(mslot, 0, 0) = trk.nInsideMinusOneHits();
    NTailMinusOneHits  (mslot, 0, 0) = trk.nTailMinusOneHits();

    LastHoT[mslot]  = trk.getLastHitOnTrack();
    CombCand[mslot] = trk.combCandidate();
  }

  void copy_out(TrackCand& trk, const int mslot, const int tslot) const
  {
    Err[tslot].CopyOut(mslot, trk.errors_nc().Array());
    Par[tslot].CopyOut(mslot, trk.parameters_nc().Array());

    trk.setCharge(Chg  (mslot, 0, 0));
    trk.setChi2  (Chi2 (mslot, 0, 0));
    trk.setLabel (Label(mslot, 0, 0));

    trk.setLastCcIndex (LastHitCcIndex(mslot, 0, 0));
    trk.setNFoundHits  (NFoundHits    (mslot, 0, 0));
    trk.setNMissingHits(NMissingHits  (mslot, 0, 0));

    trk.setNInsideMinusOneHits(NInsideMinusOneHits(mslot, 0, 0));
    trk.setNTailMinusOneHits  (NTailMinusOneHits  (mslot, 0, 0));

    trk.setCombCandidate( CombCand[mslot] );
  }

  void add_hit(const int mslot, int index, int layer)
  {
    int &n_tot_hits = NHits(mslot, 0, 0);
    int &n_fnd_hits = NFoundHits(mslot, 0, 0);

    if (n_tot_hits < Config::nMaxTrkHits)
    {
      HoTArrs[mslot][n_tot_hits++] = { index, layer };
      if (index >= 0) { ++n_fnd_hits; }
    }
    else
    {
      // printf("WARNING MkFinder::add_hit hit-on-track limit reached for label=%d\n", label_);

      const int LH = Config::nMaxTrkHits - 1;

      if (index >= 0)
      {
        if (HoTArrs[mslot][LH].index < 0)
          ++n_fnd_hits;
        HoTArrs[mslot][LH] = { index, layer };
      }
      else if (index == -2)
      {
        if (HoTArrs[mslot][LH].index >= 0)
          --n_fnd_hits;
        HoTArrs[mslot][LH] = { index, layer };
      }
    }
  }

  int num_all_minus_one_hits(const int mslot) const
  {
    return NInsideMinusOneHits(mslot, 0, 0) + NTailMinusOneHits(mslot, 0, 0);
  }

  int num_inside_minus_one_hits(const int mslot) const
  {
    return NInsideMinusOneHits(mslot, 0, 0);
  }
};

} // end namespace mkfit
#endif
