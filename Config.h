#ifndef _config_
#define _config_

#include <algorithm>

#define PRINTOUTS_FOR_PLOTS

namespace Config
{
  // math general --> from namespace TMath
  constexpr float    PI    = 3.14159265358979323846;
  constexpr float TwoPI    = 6.28318530717958647692;
  constexpr float PIOver2  = Config::PI / 2.0;
  constexpr float PIOver4  = Config::PI / 4.0;
  constexpr float InvPI    = 1.0 / Config::PI;
  constexpr float RadToDeg = 180.0 / Config::PI;
  constexpr float DegToRad = Config::PI / 180.0;
  constexpr double Sqrt2   = 1.4142135623730950488016887242097;

  // general parameters of matrices
  constexpr int nParams = 6;

  // config on main + mkFit
  extern int nTracks; //defined in Config.cc by default or when reading events from file
  extern int nEvents;

  // config on main -- for geometry
  constexpr int   nLayers   = 10;
  constexpr float fRadialSpacing   = 4.;
  constexpr float fRadialExtent    = 0.01;
  constexpr float fInnerSensorSize = 5.0; // approximate sensor size in cm
  constexpr float fOuterSensorSize = Config::fInnerSensorSize * 2.;
  constexpr float fEtaDet          = 1;  // 1 from chep

  constexpr float cmsAvgRads[10] = {4.42,7.31,10.17,25.65,33.81,41.89,49.67,60.95,69.11,78.19}; // cms average radii
  constexpr float cmsDeltaRad = 2.5; //fixme! using constant 2.5 cm, to be taken from layer properties

  // config on Event
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
  constexpr float hitposerrR  = Config::hitposerrXY / 10.;
  constexpr float varXY       = Config::hitposerrXY * Config::hitposerrXY;
  constexpr float varZ        = Config::hitposerrZ  * Config::hitposerrZ;
  constexpr float varR        = Config::hitposerrR  * Config::hitposerrR;

  constexpr int nTotHit = Config::nLayers; // for now one hit per layer for sim

  // scattering simulation
  constexpr float X0 = 9.370; // cm, from http://pdg.lbl.gov/2014/AtomicNuclearProperties/HTML/silicon_Si.html // Pb = 0.5612 cm
  constexpr float xr = 0.1; //  -assumes radial impact. This is bigger than what we have in main --> shouldn't it be the parameter below??? if radial impact??
  //const     float xr = std::sqrt(Config::beamspotX*Config::beamspotX + Config::beamspotY*Config::beamspotY); 

  // Config for seeding
  constexpr int   nlayers_per_seed = 3;
  constexpr float chi2seedcut  = 9.0;
  constexpr float alphaBeta    = 0.0520195; // 0.0458378 --> for d0 = .0025 cm --> analytically derived... depends on geometry of detector --> from mathematica
  constexpr float dEtaSeedTrip = 0.6; // for almost max efficiency --> empirically derived... depends on geometry of detector
  constexpr float dPhiSeedTrip = 0.0458712; // numerically+semianalytically derived... depends on geometry of detector

  // Config for propagation
  constexpr int Niter = 5;
  constexpr float Bfield = 3.8112;
  constexpr bool doIterative = true;
  constexpr bool useSimpleJac = false;
  constexpr bool useCurvJac   = false;
  constexpr bool useTrigApprox = true;

  // Config for Hit and BinInfoUtils
  constexpr int   nPhiPart   = 1260;
  constexpr float fPhiFactor = nPhiPart / TwoPI;
  constexpr int   nEtaPart   = 11;
  constexpr int   nEtaBin    = 2 * nEtaPart - 1;

  constexpr float        fEtaFull  = 2 * Config::fEtaDet;
  constexpr float        lEtaPart  = Config::fEtaFull/float(Config::nEtaPart);
  constexpr float        lEtaBin   = Config::lEtaPart/2.;
  constexpr float        fEtaOffB1 = Config::fEtaDet;
  constexpr float        fEtaFacB1 = Config::nEtaPart / Config::fEtaFull;
  constexpr float        fEtaOffB2 = Config::fEtaDet - Config::fEtaFull / (2 * Config::nEtaPart);
  constexpr float        fEtaFacB2 = (Config::nEtaPart>1 ? (Config::nEtaPart - 1) / (Config::fEtaFull - Config::fEtaFull / Config::nEtaPart) : 1);

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

  // matrix config
  // XXXX MT this should be renamed, made constexpr
#ifndef MAX_HITS
#define MAX_HITS 10
#endif

  //fixme: these should not be constant and modified when nTracks is set from reading a file
  constexpr int maxHitsConsidered = 25;
  const     int maxHitsPerBunch   = std::max(100, nTracks * 12 / 10 / nEtaPart) + maxHitsConsidered;

  constexpr int maxCandsPerSeed   = 6;
  constexpr int maxHolesPerCand   = 2;
  const     int maxCandsPerEtaBin = std::max(100, maxCandsPerSeed * nTracks / nEtaPart);
  // Effective eta bin is one half of nEtaPart -- so the above is twice the "average".
  // Note that last and first bin are 3/4 nEtaPart ... but can be made 1/4 by adding
  // additional bins on each end.

  // XXX std::min/max have constexpr versions in c++-14.

  extern int    numThreadsFinder;
  extern int    numThreadsSimulation;

  extern bool   clonerUseSingleThread;
  extern int    finderReportBestOutOfN;

  extern bool   useCMSGeom;
  extern bool   readCmsswSeeds;

  const std::string inputFile = "cmssw.simtracks.SingleMu1GeV.10k.new.txt";
  //const std::string inputFile = "cmssw.simtracks.SingleMu10GeV.10k.new.txt";
  //const std::string inputFile = "cmssw.rectracks.SingleMu1GeV.10k.new.txt";
  //const std::string inputFile = "cmssw.rectracks.SingleMu10GeV.10k.new.txt";

#ifdef USE_MATRIPLEX

  #ifndef MPT_SIZE
    #ifdef __MIC__
      #define MPT_SIZE 16
    #else
      #define MPT_SIZE 8
    #endif
  #endif

  #ifndef THREAD_BINDING
  #define THREAD_BINDING spread
  #endif

#endif
};

#endif 
