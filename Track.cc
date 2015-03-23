#include "Track.h"

// find the simtrack that provided the most hits
SimTkIDInfo Track::SimTrackIDInfo() const
{
  if (hits_.size() == 0){
    return SimTkIDInfo(0,0);
  }

  std::vector<unsigned int> mctrack;
  for (auto&& ihit : hits_){
    mctrack.push_back(ihit.mcTrackID());
  }
  std::sort(mctrack.begin(), mctrack.end()); // ensures all elements are checked properly

  auto mtrk(mctrack[0]), mcount(0U), m(mctrack[0]), c(0U);

  for (auto i : mctrack) {
    if (i == m) { ++c; } else { c = 0; }
    if (c >= mcount) { mtrk = m; mcount = c; }
    m = i;
  }

  return SimTkIDInfo (mtrk,mcount);
}

