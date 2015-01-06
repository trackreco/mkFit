#include "Track.h"

unsigned int getPhiPartition(float phi) {
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<TMath::Pi())) std::cout << "anomalous phi=" << phi << std::endl;
  float phiPlusPi  = phi+TMath::Pi();
  unsigned int bin = phiPlusPi*10;
  return bin;
}

unsigned int getEtaPartition(float eta, float etaDet){
  float etaPlusEtaDet  = eta + etaDet;
  float twiceEtaDet    = 2.0*etaDet;
  unsigned int bin     = (etaPlusEtaDet * 10.) / twiceEtaDet;  // ten bins for now ... update if changed in Event.cc
  return bin;
}

float getPhi(float x, float y) { 
  return std::atan2(y,x); 
}

float getEta(float x, float y, float z) {
  float theta = atan2( std::sqrt(x*x+y*y), z );
  return -1. * log( tan(theta/2.) );
}

/*unsigned int getZPartition(float z, float zPlane){
  float zPlusPlane  = z + zPlane;
  float zPartSize   = 2*zPlane / 10.; // Hardcode in 10 partitions --> need to change vector of vectors if 10 changes
  unsigned int bin  = zPlusPlane / zPartSize; 
  return bin;
}
*/

// find the simtrack that provided the most hits
std::pair<unsigned int,unsigned int> Track::SimTrackIDInfo() const
{
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

  std::pair<unsigned int, unsigned int> simIDInfo (mtrk,mcount);
  return simIDInfo;
}

