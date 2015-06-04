#ifndef _bininfoutils_
#define _bininfoutils_

#include "TMath.h"

namespace Config {
  static constexpr const float    PI = 3.14159265358979323846;
  static constexpr const float TwoPI = 6.28318530717958647692;

  static constexpr const unsigned int nPhiPart   = 63; //80; 
  static constexpr const float        nPhiFactor = nPhiPart / TwoPI;
  static constexpr const unsigned int nEtaPart   = 10; // 11;
  static constexpr const float        fEtaDet    = 2;  // 1
  static constexpr const unsigned int nEtaBin    = 2*nEtaPart - 1;

  static constexpr const float        fEtaFull  = 2 * fEtaDet;
  static constexpr const float        lEtaPart  = fEtaFull/float(nEtaPart);
  static constexpr const float        lEtaBin   = lEtaPart/2.;
  static constexpr const float        fEtaOffB1 = fEtaDet;
  static constexpr const float        fEtaFacB1 = nEtaPart / fEtaFull;
  static constexpr const float        fEtaOffB2 = fEtaDet - fEtaFull / (2 * nEtaPart);
  static constexpr const float        fEtaFacB2 = (nEtaPart - 1) / (fEtaFull - fEtaFull / nEtaPart);

  // This is for extra bins narrower ... thinking about this some more it
  // seems it would be even better to have two more exta bins, hanging off at
  // both ends.
  //
  // Anyway, it doesn't matter ... as with wide vertex region this eta binning
  // won't make much sense -- will have to be done differently for different
  // track orgin hypotheses. In about a year or so.

  inline int getEtaBin(float eta)
  {

    //in this case we are out of bounds
    if (fabs(eta)>fEtaDet) return -1;

    //first and last bin have extra width
    if (eta<(lEtaBin-fEtaDet)) return 0;
    if (eta>(fEtaDet-lEtaBin)) return nEtaBin-1;

    //now we can treat all bins as if they had same size
    return int( (eta+fEtaDet-lEtaBin/2.)/lEtaBin );

  }

  inline int getBothEtaBins(float eta, int& b1, int& b2)
  {
    b1 = b2 = -1;

    if (eta < -fEtaDet || eta > fEtaDet)
    {
      return 0;
    }

    int b1p = std::floor((eta + fEtaOffB1) * fEtaFacB1);
    int b2p = std::floor((eta + fEtaOffB2) * fEtaFacB2);

    // printf("b1' = %d   b2' = %d\n", b1p, b2p);

    int cnt = 0;
    if (b1p >= 0 && b1p < nEtaPart)
    {
      b1 = 2 * b1p;
      ++cnt;
    }
    if (b2p >= 0 && b2p < nEtaPart - 1)
    {
      b2 = 2 * b2p + 1;
      ++cnt;
    }

    // printf("b1  = %d   b2  = %d\n", b1, b2);

    return cnt;
  }
};

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
