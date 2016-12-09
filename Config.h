#ifndef _config_
#define _config_

#include <algorithm>
#include <cmath>
#include <string> // won't compile on clang gcc for mac OS w/o this!

//#define PRINTOUTS_FOR_PLOTS
#define CCSCOORD

namespace Config
{
  // super debug mode in SMatrix
  extern bool super_debug;

  // math general --> from namespace TMath
  constexpr float    PI    = 3.14159265358979323846;
  constexpr float TwoPI    = 6.28318530717958647692;
  constexpr float PIOver2  = Config::PI / 2.0f;
  constexpr float PIOver4  = Config::PI / 4.0f;
  constexpr float PI3Over4 = 3.0f * Config::PI / 4.0f;
  constexpr float InvPI    = 1.0f / Config::PI;
  constexpr float RadToDeg = 180.0f / Config::PI;
  constexpr float DegToRad = Config::PI / 180.0f;
  constexpr double Sqrt2   = 1.4142135623730950488016887242097;
  constexpr float sol      = 0.299792458; // speed of light in nm/s

  // general parameters of matrices
  constexpr int nParams = 6;

  // config on main + mkFit
  extern int nTracks; //defined in Config.cc by default or when reading events from file
  extern int nEvents;

  // config on main -- for geometry
  constexpr int   nLayers   = 10; // default: 10; cmssw tests: 13, 17, 26 (for endcap)
  constexpr float fRadialSpacing   = 4.;
  constexpr float fRadialExtent    = 0.01;
  constexpr float fInnerSensorSize = 5.0; // approximate sensor size in cm
  constexpr float fOuterSensorSize = Config::fInnerSensorSize * 2.;
  constexpr float fEtaDet          = 1; // default: 1; cmssw tests: 2, 2.5

  //constexpr float cmsAvgRads[13] = {4.42,7.31,10.17,25.58,33.98,41.79,49.78,60.78,69.2,77.96,86.80,96.53,108.00}; // cms average radii, noSplit version
  constexpr float cmsAvgRads[17] = {4.42,7.31,10.17,25.58,25.58,33.98,33.98,41.79,49.78,60.57,61.00,69.41,68.98,77.96,86.80,96.53,108.00}; // cms average radii, split version
  constexpr float cmsDeltaRad = 2.5; //fixme! using constant 2.5 cm, to be taken from layer properties

  constexpr float cmsAvgZs[26]         = {35.5,48.5,79.8,79.8,92.6,92.6,105.6,105.6,131.3,131.3,145.3,145.3,159.3,159.3,173.9,173.9,187.8,187.8,205.4,205.4,224.0,224.0,244.4,244.4,266.3,266.3}; // cms average z
  constexpr float cmsDiskMinRs[26]     = { 5.7, 5.7,23.1,22.8,23.1,22.8, 23.1, 22.8, 23.3, 23.0, 23.3, 23.0, 23.3, 23.0, 31.6, 34.4, 31.6, 34.4, 31.6, 34.4, 59.9, 38.8, 59.9, 38.8, 59.9, 49.9};
  constexpr float cmsDiskMaxRs[26]     = {14.7,14.7,50.8,42.0,50.8,42.0, 50.8, 42.0, 76.1,110.0, 76.1,110.0, 76.1,110.0, 75.9,109.7, 75.9,109.7, 75.9,109.7, 75.9,109.4, 75.9,109.4, 75.9,109.4};
  constexpr float cmsDiskMinRsHole[26] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0, 42.0,  0.0, 42.0,  0.0, 42.0,  0.0, 42.1,  0.0, 42.1,  0.0, 42.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};
  constexpr float cmsDiskMaxRsHole[26] = {999.,999.,999.,999.,999.,999., 999., 999., 59.9, 999., 59.9, 999., 59.9, 999., 59.7, 999., 59.7, 999., 59.7, 999., 999., 999., 999., 999., 999., 999.};
  const     float g_disk_dr[]        = {   1,   1,  20,  20,  20,  20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20};

  const float g_layer_zwidth[] = { 10, 14, 18, 23, 28, 32, 37, 42, 48, 52 }; //default
  const float g_layer_dz[]     = { 0.6, 0.55, 0.5, 0.5, 0.45, 0.4, 0.4, 0.4, 0.35, 0.35 }; //default

  /* const float g_layer_zwidth[] = { 30, 30, 30, 70, 70, 70, 70, 70, 70, 110, 110, 110, 110, 110, 110, 110, 110 }; //cmssw tests */
  /* const float g_layer_dz[] = { 1, 1, 1, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20 }; //cmssw tests */

  // config for material effects in cmssw
  constexpr float rangeZME = 300.;
  constexpr int   nBinsZME = 300;
  constexpr float rangeRME = 110.;
  constexpr int   nBinsRME = 110;

  constexpr float Rl[126] =
  {
    0.026,0.020,0.020,0.020,0.025,0.109,0.068,0.025,0.035,0.023,
    0.040,0.040,0.060,0.036,0.050,0.180,0.090,0.070,0.096,0.048,
    0.022,0.040,0.060,0.012,0.025,0.030,0.015,0.040,0.045,0.015,
    0.040,0.040,0.015,0.040,0.040,0.015,0.040,0.040,0.103,0.145,
    0.279,0.103,0.145,0.279,0.103,0.145,0.279,0.103,0.145,0.279,
    0.036,0.071,0.036,0.071,0.036,0.071,0.036,0.071,0.036,0.071,
    0.036,0.071,0.027,0.054,0.027,0.054,0.027,0.054,0.027,0.054,
    0.027,0.054,0.027,0.054,0.047,0.024,0.047,0.024,0.047,0.047,
    0.024,0.047,0.024,0.047,0.051,0.025,0.051,0.025,0.051,0.051,
    0.025,0.051,0.025,0.051,0.053,0.026,0.053,0.026,0.053,0.053,
    0.026,0.053,0.026,0.053,0.047,0.024,0.047,0.047,0.047,0.024,
    0.047,0.047,0.045,0.023,0.045,0.045,0.045,0.023,0.045,0.045,
    0.045,0.023,0.045,0.045,0.023,0.045
  };

  constexpr float Xi[126] =
  {
    0.054e-03,0.037e-03,0.040e-03,0.037e-03,0.040e-03,0.200e-03,0.130e-03,0.057e-03,0.080e-03,0.048e-03,
    0.060e-03,0.085e-03,0.120e-03,0.075e-03,0.110e-03,0.340e-03,0.170e-03,0.070e-03,0.220e-03,0.110e-03,
    0.056e-03,0.070e-03,0.120e-03,0.024e-03,0.040e-03,0.050e-03,0.034e-03,0.070e-03,0.070e-03,0.033e-03,
    0.070e-03,0.070e-03,0.033e-03,0.070e-03,0.070e-03,0.033e-03,0.070e-03,0.070e-03,0.270e-03,0.270e-03,
    0.480e-03,0.240e-03,0.390e-03,0.550e-03,0.270e-03,0.270e-03,0.480e-03,0.240e-03,0.390e-03,0.550e-03,
    0.080e-03,0.160e-03,0.080e-03,0.160e-03,0.080e-03,0.160e-03,0.080e-03,0.160e-03,0.080e-03,0.160e-03,
    0.080e-03,0.160e-03,0.060e-03,0.120e-03,0.060e-03,0.120e-03,0.060e-03,0.120e-03,0.060e-03,0.120e-03,
    0.060e-03,0.120e-03,0.060e-03,0.120e-03,0.110e-03,0.050e-03,0.110e-03,0.050e-03,0.110e-03,0.110e-03,
    0.050e-03,0.110e-03,0.050e-03,0.110e-03,0.120e-03,0.060e-03,0.120e-03,0.060e-03,0.120e-03,0.120e-03,
    0.060e-03,0.120e-03,0.060e-03,0.120e-03,0.120e-03,0.060e-03,0.120e-03,0.060e-03,0.120e-03,0.120e-03,
    0.060e-03,0.120e-03,0.060e-03,0.120e-03,0.110e-03,0.050e-03,0.110e-03,0.110e-03,0.110e-03,0.050e-03,
    0.110e-03,0.110e-03,0.100e-03,0.050e-03,0.100e-03,0.100e-03,0.100e-03,0.050e-03,0.100e-03,0.100e-03,
    0.050e-03,0.100e-03,0.100e-03,0.050e-03,0.100e-03,0.100e-03
  };

  extern float RlgridME[Config::nBinsZME][Config::nBinsRME];
  extern float XigridME[Config::nBinsZME][Config::nBinsRME];

  // These could be parameters, layer dependent.
  static constexpr int   m_nphi = 1024;
  static constexpr float m_max_dz   = 1; // default: 1; cmssw tests: 20
  static constexpr float m_max_dphi = 0.02; // default: 0.02; cmssw tests: 0.2

  // config on Event
  constexpr float chi2Cut = 15.;// default: 15.; cmssw tests: 30.
  constexpr float nSigma  = 3.;
  constexpr float minDPhi = 0.;// default: 0.;  cmssw tests: 0.01;
  constexpr float maxDPhi = Config::PI;
  constexpr float minDEta = 0.;
  constexpr float maxDEta = 1.0;
  constexpr float minDZ = 0.; // default: 0.; cmssw tests: 10.;

  // Configuration for simulation info
  //CMS beam spot width 25um in xy and 5cm in z 
  constexpr float beamspotX = 0.1;
  constexpr float beamspotY = 0.1;
  constexpr float beamspotZ = 1.0;
  
  constexpr float minSimPt = 0.5;
  constexpr float maxSimPt = 10.;

  constexpr float maxEta   = 1.0;

  constexpr float hitposerrXY = 0.01; // resolution is 100um in xy --> more realistic scenario is 0.003
  constexpr float hitposerrZ  = 0.1; // resolution is 1mm in z
  constexpr float hitposerrR  = Config::hitposerrXY / 10.0f;
  constexpr float varXY       = Config::hitposerrXY * Config::hitposerrXY;
  constexpr float varZ        = Config::hitposerrZ  * Config::hitposerrZ;
  constexpr float varR        = Config::hitposerrR  * Config::hitposerrR;

  constexpr int nTotHit = Config::nLayers; // for now one hit per layer for sim

  // scattering simulation
  constexpr float X0 = 9.370; // cm, from http://pdg.lbl.gov/2014/AtomicNuclearProperties/HTML/silicon_Si.html // Pb = 0.5612 cm
  constexpr float xr = 0.1; //  -assumes radial impact. This is bigger than what we have in main --> shouldn't it be the parameter below??? if radial impact??
  //const     float xr = std::sqrt(Config::beamspotX*Config::beamspotX + Config::beamspotY*Config::beamspotY); 

  // Config for seeding
  extern    int   nlayers_per_seed;
  constexpr float chi2seedcut  = 9.0;
  constexpr float lay01angdiff = 0.0634888; // analytically derived... depends on geometry of detector --> from mathematica ... d0 set to one sigma of getHypot(bsX,bsY)
  constexpr float lay02angdiff = 0.11537;
  constexpr float dEtaSeedTrip = 0.06; // for almost max efficiency --> empirically derived... depends on geometry of detector
  constexpr float dPhiSeedTrip = 0.0458712; // numerically+semianalytically derived... depends on geometry of detector
  static const float seed_z2cut= (nlayers_per_seed * fRadialSpacing) / std::tan(2.0f*std::atan(std::exp(-1.0f*dEtaSeedTrip)));
  constexpr float seed_z0cut   = beamspotZ * 3.0f; // 3cm
  constexpr float seed_z1cut   = hitposerrZ * 3.6f; // 3.6 mm --> to match efficiency from chi2cut
  constexpr float seed_d0cut   = 0.5f; // 5mm
  extern bool cf_seeding;
  extern bool findSeeds;

  // Config for propagation
  constexpr int Niter = 5;
  constexpr bool useTrigApprox = true;

  // Config for Bfield
  constexpr float Bfield = 3.8112;
  constexpr float mag_c1 = 3.8114;
  constexpr float mag_b0 = -3.94991e-06;
  constexpr float mag_b1 = 7.53701e-06;
  constexpr float mag_a  = 2.43878e-11;
  
  // Config for seeding as well... needed bfield
  constexpr float maxCurvR = (100 * minSimPt) / (sol * Bfield); // in cm

  // Config for Hit and BinInfoUtils
  constexpr int   nPhiPart   = 1260;
  constexpr float fPhiFactor = nPhiPart / TwoPI;
  constexpr int   nEtaPart   = 11;  // 1 is better for GPU best_hit
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
  constexpr float blowupfit = 10.0;
  constexpr float ptinverr049 = 0.0078; // 0.0075; // errors used for MC only fit, straight from sim tracks, outward with simple geometry
  constexpr float phierr049   = 0.0017; // 0.0017;
  constexpr float thetaerr049 = 0.0033; // 0.0031; 
  // parameters for layers 0,1,2 // --> ENDTOEND with "real seeding", fit is outward by definition, with poly geo
  constexpr float ptinverr012 = 0.12007; // 0.1789;  -->old values from only MC seeds
  constexpr float phierr012   = 1.0; // found empirically 0.00646; // 0.0071 
  constexpr float thetaerr012 = 0.2; // also found empirically 0.01366; // 0.0130; 

  // config on fitting
  extern bool cf_fitting;

  //fixme: these should not be constant and modified when nTracks is set from reading a file
  constexpr int maxHitsConsidered = 25;
  extern    int maxHitsPerBunch;

  constexpr int maxCandsPerSeed   = 6; //default: 6; cmssw tests: 3
  constexpr int maxHolesPerCand   = 2;
  extern    int maxCandsPerEtaBin;

  // config on validation
  extern bool normal_val;
  extern bool full_val;
  extern bool fit_val;

  // Effective eta bin is one half of nEtaPart -- so the above is twice the "average".
  // Note that last and first bin are 3/4 nEtaPart ... but can be made 1/4 by adding
  // additional bins on each end.

  // XXX std::min/max have constexpr versions in c++-14.
  extern int    numThreadsFinder;
  extern int    numThreadsSimulation;

  // For GPU computations
  extern int    numThreadsEvents;
  extern int    numThreadsReorg;

  extern bool   clonerUseSingleThread;
  extern int    finderReportBestOutOfN;

  extern int    numSeedsPerTask;

  // number of layer1 hits for finding seeds per task
  extern int    numHitsPerTask;
  
  extern bool   useCMSGeom;
  extern bool   readCmsswSeeds;

  extern bool   endcapTest;

  const std::string inputFile = "cmssw.simtracks.SingleMu1GeV.10k.new.txt";
  //const std::string inputFile = "cmssw.simtracks.SingleMu10GeV.10k.new.txt";
  //const std::string inputFile = "cmssw.rectracks.SingleMu1GeV.10k.new.txt";
  //const std::string inputFile = "cmssw.rectracks.SingleMu10GeV.10k.new.txt";

  void RecalculateDependentConstants();
  
  inline float BfieldFromZR(const float z, const float r)
  {
    return (Config::mag_b0*z*z + Config::mag_b1*z + Config::mag_c1)*(Config::mag_a*r*r + 1.f);
  }

#ifdef USE_MATRIPLEX

  #ifndef MPT_SIZE
    #ifdef __MIC__
      #define MPT_SIZE 16
    #elif defined USE_CUDA
      #define MPT_SIZE 8
    #else
      #define MPT_SIZE 8
    #endif
  #endif

  #ifndef THREAD_BINDING
  #define THREAD_BINDING spread
  #endif

#endif

#if defined(__CUDACC__)
  #define CUDA_CALLABLE __host__ __device__
#else
  #define CUDA_CALLABLE 
#endif
};

#endif 
