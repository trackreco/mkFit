#ifndef _bininfoutils_
#define _bininfoutils_

#include "TMath.h"

namespace Config {
  static constexpr const unsigned int nPhiPart   = 63;
  static constexpr const unsigned int nPhiFactor = 10;  // nPhiPart/2pi
  static constexpr const unsigned int nEtaPart   = 10;
  static constexpr const float        etaDet     = 2.0;
};

typedef std::pair<unsigned int,unsigned int> BinInfo;
// phi-eta partitioning map: vector of vector of vectors of std::pairs. 
// vec[nLayers][nEtaBins][nPhiBins]
typedef std::vector<std::vector<std::vector<BinInfo> > > BinInfoMap;
typedef std::vector<unsigned int> hitIndices;
typedef hitIndices::iterator hitIdxIter;

inline unsigned int getPhiPartition(float phi){
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<TMath::Pi())) std::cout << "anomalous phi=" << phi << std::endl;
  const float phiPlusPi  = phi+TMath::Pi();
  const unsigned int bin = phiPlusPi*Config::nPhiFactor;
  return bin;
}

inline unsigned int getEtaPartition(float eta){
  const float etaPlusEtaDet  = eta + Config::etaDet;
  const unsigned int bin     = (etaPlusEtaDet * Config::nEtaPart) / (2.0 * Config::etaDet);  // ten bins for now ... update if changed in Event.cc
  return bin;
}

inline float normalizedPhi(float phi) {
  return std::fmod(phi, (float) M_PI);
}

#ifdef ETASEG
inline float normalizedEta(float eta) {
  if (eta < -Config::etaDet ) eta = -Config::etaDet+.00001;
  if (eta >  Config::etaDet ) eta =  Config::etaDet-.00001;
  return eta;
}
#endif

hitIndices getCandHitIndices(const float, const float, const float, const float, const unsigned int, const BinInfoMap&);

#endif
