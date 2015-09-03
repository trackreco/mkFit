#ifndef _config_
#define _config_

#include <algorithm>

namespace Config{  
  // math general --> from namespace TMath
  static constexpr const float    PI    = 3.14159265358979323846;
  static constexpr const float TwoPI    = 6.28318530717958647692;
  static constexpr const float PIOver2  = Config::PI / 2.0;
  static constexpr const float PIOver4  = Config::PI / 4.0;
  static constexpr const float InvPI    = 1.0 / Config::PI;
  static constexpr const float RadToDeg = 180.0 / Config::PI;
  static constexpr const float DegToRad = Config::PI / 180.0;
  static constexpr const double Sqrt2   = 1.4142135623730950488016887242097;

  // config on main + mkFit
  static constexpr const unsigned int nTracks = 20000;
  static constexpr const unsigned int nEvents = 10;

  // config on main -- for geometry
  static constexpr const unsigned int nLayers   = 10;
  static constexpr const float fRadialSpacing   = 4.;
  static constexpr const float fRadialExtent    = 0.01;
  static constexpr const float fInnerSensorSize = 5.0; // approximate sensor size in cm
  static constexpr const float fOuterSensorSize = fInnerSensorSize * 2.;
  static constexpr const float fEtaDet          = 1;  // 1 from chep

  // config on Event
  static constexpr const unsigned int maxCand = 10;
  static constexpr const float chi2Cut = 15.;
  static constexpr const float nSigma = 3.;
  static constexpr const float minDPhi = 0.;
  static constexpr const float maxDPhi = Config::PI;
  static constexpr const float minDEta = 0.;
  static constexpr const float maxDEta = 1.0;

  // Configuration for simulation info
  //CMS beam spot width 25um in xy and 5cm in z 
  static constexpr const float beamspotX = 0.1;
  static constexpr const float beamspotY = 0.1;
  static constexpr const float beamspotZ = 1.0;
  
  static constexpr const float minSimPt = 0.5;
  static constexpr const float maxSimPt = 10.;

  static constexpr const float maxEta   = 1.0;

  static constexpr const float hitposerrXY = 0.01; // resolution is 100um in xy
  static constexpr const float hitposerrZ  = 0.1; // resolution is 1mm in z
  static constexpr const float hitposerrR  = hitposerrXY/10.;
  static constexpr const float varXY = hitposerrXY*hitposerrXY;
  static constexpr const float varZ  = hitposerrZ*hitposerrZ;
  static constexpr const float varR  = hitposerrR*hitposerrR;

  static constexpr const unsigned int nTotHit = nLayers; // for now one hit per layer for simu

  // Config for seeding
  static constexpr const unsigned int nlayers_per_seed = 3;

  static constexpr const float chi2seedcut = 9.0;

  // Config for Hit and BinInfoUtils
  static constexpr const unsigned int nPhiPart   = 80; //315 = 63*5  //63; //80; 
  static constexpr const float        fPhiFactor = nPhiPart / TwoPI;
  static constexpr const unsigned int nEtaPart   = 11; //11;
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

  // Config for Conformal fitter --> these change depending on inward/outward, which tracks used (MC vs reco), geometry, layers used, track params generated...
  // parameters for layers 0,4,9 
  static constexpr const float ptinverr049 = 0.0078; // 0.0075; // errors used for MC only fit, straight from sim tracks, outward with simple geometry
  static constexpr const float phierr049   = 0.0017; // 0.0017;
  static constexpr const float thetaerr049 = 0.0033; // 0.0031; 
  // parameters for layers 0,1,2 // --> ENDTOEND with "real seeding", fit is outward by definition, with poly geo
  static constexpr const float ptinverr012 = 0.1234; // 0.1789;  -->old values from only MC seeds
  static constexpr const float phierr012   = 0.0071; // 0170; 
  static constexpr const float thetaerr012 = 0.0130; // 0.0137; 

  // config for MPlex HitStructures 
  static const unsigned int g_MaxHitsPerBunch = std::max( (unsigned int)(100), (unsigned int)((Config::nTracks * 2) / Config::nEtaPart));

  static constexpr const unsigned int g_MaxCandsPerSeed   = 6;
  static const unsigned int g_MaxCandsPerEtaBin = std::max( (unsigned int)(100), (unsigned int)((Config::g_MaxCandsPerSeed * Config::nTracks) / Config::nEtaPart));
  // Effective eta bin is one half of nEtaPart -- so the above is twice the "average".
  // Note that last and first bin are 3/4 nEtaPart ... but can be made 1/4 by adding
  // additional bins on each end.

  //matrix config
#ifndef MAX_HITS
#define MAX_HITS 10
#endif

#ifdef USE_MATRIPLEX

  #ifndef MPT_SIZE
    #ifdef __MIC__
      #define MPT_SIZE 16
    #else
      #define MPT_SIZE 8
    #endif
  #endif

  #ifndef NUM_THREADS
  #define NUM_THREADS 1
  #endif

  #ifndef NUM_THREADS_SIM
  #define NUM_THREADS_SIM 60
  #endif

  #ifndef THREAD_BINDING
  #define THREAD_BINDING spread
  #endif

#endif
};

#endif 
