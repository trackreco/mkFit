#ifndef _config_
#define _config_

#include <algorithm>

namespace Config{  
  // math general --> from namespace TMath
  constexpr float    PI    = 3.14159265358979323846;
  constexpr float TwoPI    = 6.28318530717958647692;
  constexpr float PIOver2  = Config::PI / 2.0;
  constexpr float PIOver4  = Config::PI / 4.0;
  constexpr float InvPI    = 1.0 / Config::PI;
  constexpr float RadToDeg = 180.0 / Config::PI;
  constexpr float DegToRad = Config::PI / 180.0;
  constexpr double Sqrt2   = 1.4142135623730950488016887242097;

  // config on main + mkFit
  constexpr unsigned int nTracks = 20000;
  constexpr unsigned int nEvents = 10;

  // config on main -- for geometry
  constexpr unsigned int nLayers   = 10;
  constexpr float fRadialSpacing   = 4.;
  constexpr float fRadialExtent    = 0.01;
  constexpr float fInnerSensorSize = 5.0; // approximate sensor size in cm
  constexpr float fOuterSensorSize = Config::fInnerSensorSize * 2.;
  constexpr float fEtaDet          = 1;  // 1 from chep

  // config on Event
  constexpr unsigned int maxCand = 10;
  constexpr float chi2Cut = 15.;
  constexpr float nSigma  = 3.;
  constexpr float minDPhi = 0.;
  constexpr float maxDPhi = Config::PI;
  constexpr float minDEta = 0.;
  constexpr float maxDEta = 1.0;

  // Configuration for simulation info
  //CMS beam spot width 25um in xy and 5cm in z 
  constexpr float beamspotX = 0.1;
  constexpr float beamspotY = 0.1;
  constexpr float beamspotZ = 1.0;
  
  constexpr float minSimPt = 0.5;
  constexpr float maxSimPt = 10.;

  constexpr float maxEta   = 1.0;

  constexpr float hitposerrXY = 0.01; // resolution is 100um in xy
  constexpr float hitposerrZ  = 0.1; // resolution is 1mm in z
  constexpr float hitposerrR  = Config::hitposerrXY/10.;
  constexpr float varXY       = Config::hitposerrXY*Config::hitposerrXY;
  constexpr float varZ        = Config::hitposerrZ*Config::hitposerrZ;
  constexpr float varR        = Config::hitposerrR*Config::hitposerrR;

  constexpr unsigned int nTotHit = Config::nLayers; // for now one hit per layer for sim

  // scattering simulation
  constexpr float X0 = 9.370; // cm, from http://pdg.lbl.gov/2014/AtomicNuclearProperties/HTML/silicon_Si.html // Pb = 0.5612 cm
  constexpr float xr = 0.1; //  -assumes radial impact. This is bigger than what we have in main --> shouldn't it be the parameter below??? if radial impact??
  //const     float xr = std::sqrt(Config::beamspotX*Config::beamspotX + Config::beamspotY*Config::beamspotY); 

  // Config for seeding
  constexpr unsigned int nlayers_per_seed = 3;
  constexpr float chi2seedcut  = 9.0;
  constexpr float alphaBeta    = 0.0520195; // 0.0458378 --> for d0 = .0025 cm --> analytically derived... depends on geometry of detector --> from mathematica
  constexpr float dEtaSeedTrip = 0.6; // for almost max efficiency --> empirically derived... depends on geometry of detector
  constexpr float dPhiSeedTrip = 0.0458712; // numerically+semianalytically derived... depends on geometry of detector

  // Config for Hit and BinInfoUtils
  constexpr unsigned int nPhiPart   = 80; //315 = 63*5  //63; //80; 
  constexpr float        fPhiFactor = Config::nPhiPart / Config::TwoPI;
  constexpr unsigned int nEtaPart   = 11; //11;
  constexpr unsigned int nEtaBin    = 2*Config::nEtaPart - 1;

  constexpr float        fEtaFull  = 2 * Config::fEtaDet;
  constexpr float        lEtaPart  = Config::fEtaFull/float(Config::nEtaPart);
  constexpr float        lEtaBin   = Config::lEtaPart/2.;
  constexpr float        fEtaOffB1 = Config::fEtaDet;
  constexpr float        fEtaFacB1 = Config::nEtaPart / Config::fEtaFull;
  constexpr float        fEtaOffB2 = Config::fEtaDet - Config::fEtaFull / (2 * Config::nEtaPart);
  constexpr float        fEtaFacB2 = (Config::nEtaPart - 1) / (Config::fEtaFull - Config::fEtaFull / Config::nEtaPart);

  // This is for extra bins narrower ... thinking about this some more it
  // seems it would be even better to have two more exta bins, hanging off at
  // both ends.
  //
  // Anyway, it doesn't matter ... as with wide vertex region this eta binning
  // won't make much sense -- will have to be done differently for different
  // track orgin hypotheses. In about a year or so.

  // Config for Conformal fitter --> these change depending on inward/outward, which tracks used (MC vs reco), geometry, layers used, track params generated...
  // parameters for layers 0,4,9 
  constexpr float ptinverr049 = 0.0078; // 0.0075; // errors used for MC only fit, straight from sim tracks, outward with simple geometry
  constexpr float phierr049   = 0.0017; // 0.0017;
  constexpr float thetaerr049 = 0.0033; // 0.0031; 
  // parameters for layers 0,1,2 // --> ENDTOEND with "real seeding", fit is outward by definition, with poly geo
  constexpr float ptinverr012 = 0.1234; // 0.1789;  -->old values from only MC seeds
  constexpr float phierr012   = 0.0071; // 0170; 
  constexpr float thetaerr012 = 0.0130; // 0.0137; 

  // config for MPlex HitStructures 
  const unsigned int g_MaxHitsPerBunch = std::max( (unsigned int)(100), (unsigned int)((Config::nTracks * 2) / Config::nEtaPart));

  constexpr unsigned int g_MaxCandsPerSeed   = 6;
  const unsigned int g_MaxCandsPerEtaBin = std::max( (unsigned int)(100), (unsigned int)((Config::g_MaxCandsPerSeed * Config::nTracks) / Config::nEtaPart));
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
