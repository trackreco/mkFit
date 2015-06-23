#ifndef _bininfoutils_
#define _bininfoutils_

#include "TMath.h"
#include "Config.h"

typedef std::pair<unsigned int,unsigned int> BinInfo;
typedef std::vector<std::vector<BinInfo> > BinInfoLayerMap;
typedef std::vector<std::vector<std::vector<BinInfo> > > BinInfoMap;

inline unsigned int getPhiPartition(float phi){
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<Config::PI)) std::cout << "anomalous phi=" << phi << std::endl;
  const float phiPlusPi  = phi+Config::PI;
  const unsigned int bin = phiPlusPi*Config::nPhiFactor;
  return bin;
}

inline unsigned int getEtaPartition(float eta){
  const float etaPlusEtaDet  = eta + Config::fEtaDet;
  const unsigned int bin     = (etaPlusEtaDet * Config::nEtaPart) / (2.0 * Config::fEtaDet);  // ten bins for now ... update if changed in Event.cc
  return bin;
}

inline float normalizedPhi(float phi) {
  return std::fmod(phi, (float) Config::PI);
}

#ifdef ETASEG
inline float normalizedEta(float eta) {
  if (fabs(eta)>Config::fEtaDet) eta = (eta>0 ? Config::fEtaDet*0.99 : -Config::fEtaDet*0.99);
  return eta;
}
#endif

std::vector<unsigned int> getCandHitIndices(const unsigned int &, const unsigned int &, const unsigned int &, const unsigned int &, const BinInfoLayerMap &);

#endif
