#include "Track.h"
#include "Debug.h"

// find the simtrack that provided the most hits
void TrackExtra::setMCTrackIDInfo(const Track& trk, const std::vector<HitVec>& layerHits, const MCHitInfoVec& globalHitInfo)
{
  std::vector<unsigned int> mctrack;
  auto hitIdx = trk.nTotalHits();
  for (int ihit = 0; ihit < hitIdx; ++ihit){
    auto mchitid = layerHits[ihit][trk.getHitIdx(ihit)].mcHitID();
    mctrack.push_back(globalHitInfo[mchitid].mcTrackID_);
  }
  std::sort(mctrack.begin(), mctrack.end()); // ensures all elements are checked properly

  auto mtrk(mctrack[0]), mcount(0U), m(mctrack[0]), c(0U);

  for (auto i : mctrack) {
    if (i == m) { ++c; } else { c = 1; }
    if (c >= mcount) { mtrk = m; mcount = c; }
    m = i;
  }
  if (2*mcount > trk.nFoundHits()){ // if more, matched track --> set id info
    mcTrackID_    = mtrk;
    nHitsMatched_ = mcount;
  } else { // fake track, id = 999 999
    mcTrackID_    = 999999;
    nHitsMatched_ = mcount;
  }
  dprint("Track " << trk.label() << " best mc track " << mtrk << " count " << mcount << "/" << trk.nFoundHits());
  // need to include protection against zero size tracks --> need ID other than 999999
}

/*
void Track::write_out(FILE *fp)
{
#if 0
  Track t = clone_for_io();
  fwrite(&t, sizeof(Track), 1, fp);

  unsigned int nh = nHits();
  fwrite(&nh, sizeof(unsigned int), 1, fp);

  fwrite(&hits_[0], sizeof(Hit), nh, fp);
#endif
}

void Track::read_in(FILE *fp)
{
#if 0
  fread(this, sizeof(Track), 1, fp);

  unsigned int nh = nHits();
  fread(&nh, sizeof(unsigned int), 1, fp);

  hits_.resize(nh);
  fread(&hits_[0], sizeof(Hit), nh, fp);
#endif
}
*/
