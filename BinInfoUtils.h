#ifndef _bininfoutils_
#define _bininfoutils_

// could #include "TMath.h", but then this is problematic for non root running.  This compiles and runs without tmath.h

#include <cmath>
#include <vector>
#include "Config.h"
#include <iostream>

typedef std::pair<unsigned int,unsigned int> BinInfo;
typedef std::vector<std::vector<BinInfo> > BinInfoLayerMap;
typedef std::vector<std::vector<std::vector<BinInfo> > > BinInfoMap;

inline float downPhi(float phi){
  while (phi >= Config::PI) {phi-=Config::TwoPI;}
  return phi;
}
	
inline float upPhi(float phi){
  while (phi <= -Config::PI) {phi+=Config::TwoPI;}
  return phi;
}

inline float normalizedPhi(float phi) {
  //  return std::fmod(phi, (float) Config::PI); // return phi +pi out of phase for |phi| beyond boundary! 
  if (std::abs(phi)>=Config::PI) {phi = (phi>0 ? downPhi(phi) : upPhi(phi));}
  return phi;
}

inline unsigned int getPhiPartition(float phi){
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<Config::PI)) std::cout << "anomalous phi=" << phi << std::endl;
  //  const float phiPlusPi  = std::fmod(phi+Config::PI,Config::TwoPI); // normaliztion done here
  const float phiPlusPi = phi+Config::PI; 
  int bin = phiPlusPi*Config::fPhiFactor;
  
  // theoretically these checks below should be taken care of by normalizedPhi, however...
  // these condition checks appeared in very bizarre corner case where propagated phi == pi != Config::PI in check of normalizedPhi (but not unexpected... comparing float point numbers)
  // i.e. delta on floating point check smaller than comparison... making what should be bin = nPhiPart - 1 instead bin = nPhiPart (out of bounds!!) ...or worse if unsigned bin < 0, bin == unsigned int max!
  if (bin<0)                      bin = 0;
  else if (bin>=Config::nPhiPart) bin = Config::nPhiPart - 1;

  return (unsigned int) bin;
}

inline unsigned int getEtaPartition(float eta){
  const float etaPlusEtaDet  = eta + Config::fEtaDet;
  int bin = (etaPlusEtaDet * Config::nEtaPart) / (Config::fEtaFull);  // ten bins for now ... update if changed in Event.cc
  // use this check instead of normalizedEta()... directly bin into first/last bins eta of range of detector
  if      (bin<0)                 bin = 0; 
  else if (bin>=Config::nEtaPart) bin = Config::nEtaPart - 1;

  return (unsigned int) bin;
}

std::vector<unsigned int> getCandHitIndices(const unsigned int &, const unsigned int &, const unsigned int &, const unsigned int &, const BinInfoLayerMap &);

#endif
