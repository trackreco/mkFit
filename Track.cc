#include "Track.h"
#include "Debug.h"

// find the simtrack that provided the most hits
void TrackExtra::setMCTrackIDInfo(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo)
{
  std::vector<int> mctrack;
  auto hitIdx = trk.nTotalHits();
  for (int ihit = 0; ihit < hitIdx; ++ihit){
    if (trk.getHitIdx(ihit) >= 0) {
      auto mchitid = layerHits[ihit][trk.getHitIdx(ihit)].mcHitID();
      mctrack.push_back(globalHitInfo[mchitid].mcTrackID());
    }
  }
  std::sort(mctrack.begin(), mctrack.end()); // ensures all elements are checked properly

  auto mtrk(mctrack[0]), mcount(0), m(mctrack[0]), c(0);

  for (auto i : mctrack) {
    if (i == m) { ++c; } else { c = 1; }
    if (c >= mcount) { mtrk = m; mcount = c; }
    m = i;
  }
  if (4*mcount > 3*trk.nFoundHits()){ // if more, matched track --> set id info
    mcTrackID_    = mtrk;
    nHitsMatched_ = mcount;
  } else { // fake track, id = 999 999
    mcTrackID_    = 999999;
    nHitsMatched_ = mcount;
  }
  dprint("Track " << trk.label() << " best mc track " << mtrk << " count " << mcount << "/" << trk.nFoundHits());
  // need to include protection against zero size tracks --> need ID other than 999999
}
