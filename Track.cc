#include "Track.h"

unsigned int getPhiPartition(float phi) {
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<TMath::Pi())) std::cout << "anomalous phi=" << phi << std::endl;
  float phiPlusPi  = phi+TMath::Pi();
  unsigned int bin = phiPlusPi*10;
  return bin;
}

std::vector<unsigned int> Track::SimTrackIDs() const
{
  std::vector<unsigned int> mctrack;
  for (auto&& ihit : hits_){
    mctrack.push_back(ihit.mcIndex());
  }
  return mctrack;
}

// find the simtrack that provided the most hits
unsigned int Track::SimTrackID() const
{
  std::vector<unsigned int> mctrack(SimTrackIDs());
  std::sort(mctrack.begin(), mctrack.end());

  auto mtrk(mctrack[0]), mcount(0U), m(mctrack[0]), c(0U);

  for (auto&& i : mctrack) {
    if (i == m) { ++c; } else { c = 0; }
    if (c >= mcount) { mtrk = m; mcount = c; }
    m = i;
  }
  return mtrk;
}
