#ifndef MkFinder_h
#define MkFinder_h

#include "MkBase.h"

#include "Track.h"

class CandCloner;
class LayerOfHits;
class SteeringParams;

//#include "Event.h"

//#include "HitStructures.h"

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


class MkFinder : public MkBase
{
public:

  static constexpr int MPlexHitIdxMax = 16;

  using MPlexHitIdx = Matriplex::Matriplex<int, MPlexHitIdxMax, 1, NN>;
  using MPlexQHoT   = Matriplex::Matriplex<HitOnTrack, 1, 1, NN>;

  struct IdxChi2List
  {
    int   trkIdx; // candidate index
    int   hitIdx; // hit index
    int   nhits;  // number of hits (used for sorting)
    float chi2;   // total chi2 (used for sorting)
  };

  //----------------------------------------------------------------------------

  MPlexQF    Chi2;
  MPlexQI    Label;   // seed index in global seed vector (for MC truth match)

  MPlexQI    NHits;
  MPlexQI    NFoundHits;
  HitOnTrack HoTArrs[NN][Config::nMaxTrkHits];

  MPlexQI    SeedIdx; // seed index in local thread (for bookkeeping at thread level)
  MPlexQI    CandIdx; // candidate index for the given seed (for bookkeeping of clone engine)

  // Hit indices into LayerOfHits to explore.
  MPlexQI     XHitSize;
  MPlexHitIdx XHitArr;

  // Hit error / parameters for hit matching, update.
  MPlexHS    msErr;
  MPlexHV    msPar;

  //============================================================================

  MkFinder() {}

  //----------------------------------------------------------------------------

  void InputTracksAndHitIdx(const std::vector<Track>& tracks,
                            int beg, int end, bool inputProp);

  void InputTracksAndHitIdx(const std::vector<Track>& tracks,
                            const std::vector<int>  &   idxs,
                            int beg, int end, bool inputProp, int mp_offset);

  void InputTracksAndHitIdx(const std::vector<std::vector<Track>>& tracks,
                            const std::vector<std::pair<int,int>>& idxs,
                            int beg, int end, bool inputProp);

  void InputTracksAndHitIdx(const std::vector<std::vector<Track>>& tracks,
                            const std::vector<std::pair<int,IdxChi2List>>& idxs,
                            int beg, int end, bool inputProp);

  void OutputTracksAndHitIdx(std::vector<Track>& tracks,
                             int beg, int end, bool outputProp) const;

  void OutputTracksAndHitIdx(std::vector<Track>& tracks,
                             const std::vector<int>& idxs,
                             int beg, int end, bool outputProp) const;

  //----------------------------------------------------------------------------

  void SelectHitIndices(const LayerOfHits &layer_of_hits, const int N_proc, bool dump=false);

  void AddBestHit(const LayerOfHits    &layer_of_hits, const int N_proc,
                  const SteeringParams &st_par);

  //----------------------------------------------------------------------------

  void FindCandidates(const LayerOfHits &layer_of_hits,
                      std::vector<std::vector<Track>>& tmp_candidates,
		      const int offset, const int N_proc,
                      const SteeringParams &st_par);

  //----------------------------------------------------------------------------

  void FindCandidatesMinimizeCopy(const LayerOfHits &layer_of_hits, CandCloner& cloner,
                                  const int offset, const int N_proc,
                                  const SteeringParams &st_par);
  void FindCandidatesMinimizeCopyEndcap(const LayerOfHits &layer_of_hits, CandCloner& cloner,
                                        const int offset, const int N_proc);

  void UpdateWithLastHit(const LayerOfHits &layer_of_hits, int N_proc);
  void UpdateWithLastHitEndcap(const LayerOfHits &layer_of_hits, int N_proc);

  void CopyOutParErr(std::vector<std::vector<Track> >& seed_cand_vec,
                     int N_proc, bool outputProp) const;

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

  void add_hit(const int mslot, int index, int layer)
  {
    HoTArrs[mslot][NHits(mslot, 0, 0)] = { index, layer };
    ++NHits(mslot, 0, 0);
    if (index >= 0) ++NFoundHits(mslot, 0, 0);
  }
};

#endif
