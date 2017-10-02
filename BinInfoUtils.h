#ifndef _bininfoutils_
#define _bininfoutils_

// XXXXMT4K This should all be removed ...
// getCandHitIndices seems to be still used.

// could #include "TMath.h", but then this is problematic for non root running.  This compiles and runs without tmath.h

#include <cmath>
#include <vector>
#include "Config.h"
#include <iostream>

typedef std::pair<int, int>                 BinInfo;
typedef std::vector<std::vector<BinInfo>>   BinInfoLayerMap;
typedef std::vector<BinInfoLayerMap>        BinInfoMap;

CUDA_CALLABLE
inline float downPhi(float phi)
{
  while (phi >= Config::PI) {phi-=Config::TwoPI;}
  return phi;
}
	
CUDA_CALLABLE
inline float upPhi(float phi)
{
  while (phi <= -Config::PI) {phi+=Config::TwoPI;}
  return phi;
}

CUDA_CALLABLE
inline float normalizedPhi(float phi)
{
  //  return std::fmod(phi, (float) Config::PI); // return phi +pi out of phase for |phi| beyond boundary! 
  if (std::abs(phi)>=Config::PI) {phi = (phi>0 ? downPhi(phi) : upPhi(phi));}
  return phi;
}

CUDA_CALLABLE
inline int getPhiPartition(float phi)
{
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<Config::PI)) std::cout << "anomalous phi=" << phi << std::endl;
  //  const float phiPlusPi  = std::fmod(phi+Config::PI,Config::TwoPI); // normaliztion done here
  const float phiPlusPi = phi+Config::PI; 
  int bin = phiPlusPi*Config::fPhiFactor;
  
  // theoretically these checks below should be taken care of by normalizedPhi, however...
  // these condition checks appeared in very bizarre corner case where propagated phi == pi != Config::PI in check of normalizedPhi (but not unexpected... comparing float point numbers)
  // i.e. delta on floating point check smaller than comparison... making what should be bin = nPhiPart - 1 instead bin = nPhiPart (out of bounds!!) ...or worse if unsigned bin < 0, bin == int max!
  if (bin<0)                      bin = 0;
  else if (bin>=Config::nPhiPart) bin = Config::nPhiPart - 1;

  return bin;
}

inline int getEtaPartition(float eta)
{
  const float etaPlusEtaDet  = eta + Config::fEtaDet;
  int bin = (etaPlusEtaDet * Config::nEtaPart) / (Config::fEtaFull);  // ten bins for now ... update if changed in Event.cc
  // use this check instead of normalizedEta()... directly bin into first/last bins eta of range of detector
  if      (bin<0)                 bin = 0; 
  else if (bin>=Config::nEtaPart) bin = Config::nEtaPart - 1;

  return bin;
}

std::vector<int> getCandHitIndices(const int &, const int &, const int &, const int &, const BinInfoLayerMap &);

#endif
