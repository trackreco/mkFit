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
  if (fabs(phi)>=Config::PI) {phi = (phi>0 ? downPhi(phi) : upPhi(phi));}
  return phi;
}

inline unsigned int getPhiPartition(float phi){
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<Config::PI)) std::cout << "anomalous phi=" << phi << std::endl;
  //  const float phiPlusPi  = std::fmod(phi+Config::PI,Config::TwoPI); // normaliztion done here
  const float phiPlusPi  = phi+Config::PI; 
  const unsigned int bin = phiPlusPi*Config::fPhiFactor;
  return bin;
}

inline unsigned int getEtaPartition(float eta){
  const float etaPlusEtaDet  = eta + Config::fEtaDet;
  const unsigned int bin     = (etaPlusEtaDet * Config::nEtaPart) / (Config::fEtaFull);  // ten bins for now ... update if changed in Event.cc
  return bin;
}
			      
#ifdef ETASEG
inline float normalizedEta(float eta) {
  if (fabs(eta)>Config::fEtaDet) {eta = (eta>0 ? Config::fEtaDet*0.99 : -Config::fEtaDet*0.99);}
  return eta;
}
#endif

std::vector<unsigned int> getCandHitIndices(const unsigned int &, const unsigned int &, const unsigned int &, const unsigned int &, const BinInfoLayerMap &);

#endif
