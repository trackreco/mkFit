#include "Track.h"

// find the simtrack that provided the most hits
void Track::setMCTrackIDInfo() 
{
  std::vector<unsigned int> mctrack;
  for (auto&& ihit : hits_){
    mctrack.push_back(ihit.mcTrackID());
  }
  std::sort(mctrack.begin(), mctrack.end()); // ensures all elements are checked properly

  auto mtrk(mctrack[0]), mcount(0U), m(mctrack[0]), c(0U);

  for (auto i : mctrack) {
    if (i == m) { ++c; } else { c = 1; }
    if (c >= mcount) { mtrk = m; mcount = c; }
    m = i;
  }

  if (4*mcount >= 3*hits_.size()){ // if more, matched track --> set id info
    mcTrackID_    = mtrk;
    nHitsMatched_ = mcount;
  }
  else{ // fake track, id = 999 999
    mcTrackID_    = 999999;
    nHitsMatched_ = mcount;
  }

  // need to include protection against zero size tracks --> need ID other than 999999
}

