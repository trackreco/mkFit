#ifndef _config_
#define _config_

// Cram this in here for now ...
class TrackerInfo;

#include <algorithm>
#include <cmath>
#include <string> // won't compile on clang gcc for mac OS w/o this!
#include <map>

#if defined(__CUDACC__)
  #define CUDA_CALLABLE __host__ __device__
#else
  #define CUDA_CALLABLE 
#endif


//------------------------------------------------------------------------------

enum PropagationFlagsEnum
{
  PF_none              = 0,

  PF_use_param_b_field = 0x1,
  PF_apply_material    = 0x2
};

struct PropagationFlags
{
  union
  {
    struct
    {
      bool use_param_b_field  : 1;
      bool apply_material     : 1;
      // Could add: bool use_trig_approx  -- now Config::useTrigApprox = true
      // Could add: int  n_iter : 8       -- now Config::Niter = 5
    };

    unsigned int _raw_;

  };

  PropagationFlags() : _raw_(0) {}

  PropagationFlags(int pfe) :
    use_param_b_field       ( pfe & PF_use_param_b_field),
    apply_material          ( pfe & PF_apply_material)
  {}
};

//------------------------------------------------------------------------------

// Enum for input seed options
enum seedOpts {simSeeds, cmsswSeeds, findSeeds};
typedef std::map<std::string,seedOpts> seedOptsMap;

// Enum for seed cleaning options
enum cleanOpts {noCleaning, cleanSeedsN2, cleanSeedsPure, cleanSeedsBadLabel};
typedef std::map<std::string,cleanOpts> cleanOptsMap;

// Enum for cmssw matching options
enum matchOpts {trkParamBased, hitBased, labelBased};
typedef std::map<std::string, matchOpts> matchOptsMap;

//------------------------------------------------------------------------------

namespace Config
{
  extern TrackerInfo TrkInfo;

  // default file version
  constexpr int FileVersion = 1;

  // math general --> from namespace TMath
  constexpr float    PI    = 3.14159265358979323846;
  constexpr float TwoPI    = 6.28318530717958647692;
  constexpr float PIOver2  = Config::PI / 2.0f;
  constexpr float PIOver4  = Config::PI / 4.0f;
  constexpr float PI3Over4 = 3.0f * Config::PI / 4.0f;
  constexpr float InvPI    = 1.0f / Config::PI;
  constexpr float RadToDeg = 180.0f / Config::PI;
  constexpr float DegToRad = Config::PI / 180.0f;
  constexpr float sol      = 0.299792458; // speed of light in nm/s
  constexpr double Sqrt2   = 1.4142135623730950488016887242097;
  constexpr double OOSqrt2 = 1.0 / Config::Sqrt2;

  // general parameters of matrices
  constexpr int nParams = 6;

  // config on main + mkFit
  extern int nTracks; //defined in Config.cc by default or when reading events from file
  extern int nEvents;
  // XXXXMT: nTracks should be thrown out ... SMatrix and Event allocate some arrays on this
  // which can be wrong for real data or in multi-event environment

  extern std::string geomPlugin;

  // config on main -- for geometry; XXXXMT to be thrown out, too
  constexpr int   nLayers   = 10; // default: 10; cmssw tests: 13, 17, 26 (for endcap)

  // New layer constants for common barrel / endcap. I'd prefer those to go
  // into some geometry definition "plugin" -- they belong more into some Geom
  // namespace, too.
  // XXXX This needs to be generalized for other geometries !
  // TrackerInfo more or less has all this information (or could have it).
  extern    int   nTotalLayers;        // To be set by geometry plugin.
  constexpr int   nMaxSimHits    = 32; // Assuming dual hit on every barrel / endcap edge -> used in tkNtuple
  constexpr int   nMaxTrkHits    = nMaxSimHits; // Used for array sizes in MkFitter and Track

  constexpr float fRadialSpacing   = 4.;
  constexpr float fRadialExtent    = 0.01;
  constexpr float fInnerSensorSize = 5.0; // approximate sensor size in cm
  constexpr float fOuterSensorSize = Config::fInnerSensorSize * 2.;
  constexpr float fEtaDet          = 1; // default: 1; cmssw tests: 2, 2.5

  constexpr float cmsDeltaRad = 2.5; //fixme! using constant 2.5 cm, to be taken from layer properties

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

  // This will become layer dependent (in bits). To be consistent with min_dphi.
  static constexpr int m_nphi = 1024;

  // config on Event
  extern    float chi2Cut; // default: 15; cmssw: 30 (set in TrackerInfo plugin)
  // the following are only used in SMatrix version
  constexpr float nSigma  = 3.;
  constexpr float minDPhi = 0.01;// default: 0.;  cmssw tests: 0.01;
  constexpr float maxDPhi = Config::PI;
  constexpr float minDEta = 0.;
  constexpr float maxDEta = 1.0;

  // Configuration for simulation info
  // CMS beam spot width 25um in xy and 5cm in z
  constexpr float beamspotX = 0.1;
  constexpr float beamspotY = 0.1;
  constexpr float beamspotZ = 1.0;

  // XXMT4K minPt was 0.5. Figure out what is the new limit for 90cm or be
  // more flexible about finding fewer hits. Or postprocess looper candidates.
  constexpr float minSimPt = 1;
  constexpr float maxSimPt = 10.;

  // XXMT Hardhack -- transition region excluded in Simulation::setupTrackByToyMC()
  constexpr float minSimEta = -2.4;
  constexpr float maxSimEta =  2.4;
  // For testing separate EC-/BRL/EC+; -2.3--1.5 / -0.9-0.9 / 1.5-2.3
  //constexpr float minSimEta =  -0.9;
  //constexpr float maxSimEta =   0.9;

  constexpr float hitposerrXY = 0.01; // resolution is 100um in xy --> more realistic scenario is 0.003
  constexpr float hitposerrZ  = 0.1;  // resolution is 1mm in z
  constexpr float hitposerrR  = Config::hitposerrXY / 10.0f; // XXMT4K ??? I don't get this ...
  constexpr float varXY       = Config::hitposerrXY * Config::hitposerrXY;
  constexpr float varZ        = Config::hitposerrZ  * Config::hitposerrZ;
  constexpr float varR        = Config::hitposerrR  * Config::hitposerrR;

  // XXMT4K OK ... what do we do with this guy? MaxTotHit / AvgTotHit ... ?
  // For now setting it to nMaxTrkLayers which is too big ... but it seems to be
  // only used for vector::reserve() ...
  constexpr int nTotHit = Config::nMaxSimHits; // for now one hit per layer for sim

  // scattering simulation
  constexpr float X0 = 9.370; // cm, from http://pdg.lbl.gov/2014/AtomicNuclearProperties/HTML/silicon_Si.html // Pb = 0.5612 cm
  constexpr float xr = 0.1; //  -assumes radial impact. This is bigger than what we have in main --> shouldn't it be the parameter below??? if radial impact??
  //const     float xr = std::sqrt(Config::beamspotX*Config::beamspotX + Config::beamspotY*Config::beamspotY); 

  // Config for seeding
  extern    int   nlayers_per_seed;         // default: 3, cms sets from geom plugin
  constexpr int   nlayers_per_seed_max = 4; // Needed for allocation of arrays on stack.
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

  // Config for propagation - could/should enter into PropagationFlags?!
  constexpr int Niter    =  5;
  constexpr int NiterSim = 10; // Can make more steps due to near volume misses.
  constexpr bool useTrigApprox = true;

  // PropagationFlags as used during finding and fitting. Defined for each Geom in its plugin.
  extern bool             finding_requires_propagation_to_hit_pos;
  extern PropagationFlags finding_inter_layer_pflags;
  extern PropagationFlags finding_intra_layer_pflags;
  extern PropagationFlags backward_fit_pflags;
  extern PropagationFlags forward_fit_pflags;
  extern PropagationFlags seed_fit_pflags;
  extern PropagationFlags pca_prop_pflags;

  // Config for Bfield. Note: for now the same for CMS-2017 and CylCowWLids.
  constexpr float Bfield = 3.8112;
  constexpr float mag_c1 = 3.8114;
  constexpr float mag_b0 = -3.94991e-06;
  constexpr float mag_b1 = 7.53701e-06;
  constexpr float mag_a  = 2.43878e-11;

  // Config for SelectHitIndices
  // Use extra arrays to store phi and q of hits.
  // MT: This would in principle allow fast selection of good hits, if
  // we had good error estimates and reasonable *minimal* phi and q windows.
  // Speed-wise, those arrays (filling AND access, about half each) cost 1.5%
  // and could help us reduce the number of hits we need to process with bigger
  // potential gains.
  extern bool usePhiQArrays;

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

  extern    int maxCandsPerSeed; // default: 6; cms: 6  (GC had 3)
  extern    int maxHolesPerCand; // default: 2; cms  12 (should be reduced)
  extern    int maxCandsPerEtaBin;

  // Selection of simtracks from CMSSW. Used in Event::clean_cms_simtracks() and MkBuilder::prep_cmsswtracks()
  constexpr int   cmsSelMinLayers = 8;
  constexpr float cmsSelMinPt     = 0.5;

  // config on validation
  constexpr int nMinFoundHits = 7;
  constexpr float minCMSSWMatchChi2[6] = {100,100,50,50,30,20};
  constexpr float minCMSSWMatchdPhi[6] = {0.2,0.2,0.1,0.05,0.01,0.005};
  constexpr int   nCMSSWMatchHitsAfterSeed = 5;
  extern bool quality_val;
  extern bool root_val;
  extern bool cmssw_val;
  extern bool fit_val;
  extern bool readSimTrackStates; // need this to fill pulls
  extern bool inclusiveShorts;
  extern bool applyCMSSWHitMatch;
  extern matchOpts cmsswMatching;

  // config on seed cleaning
  constexpr int minNHits_seedclean = 4;
  constexpr float track1GeVradius = 87.6; // = 1/(c*B)
  constexpr float c_etamax_brl = 0.9;
  constexpr float c_dpt_brl_0 = 0.025;
  constexpr float c_dpt_ec_0 = 0.0125;
  constexpr float c_ptmax_0 = 2.0;
  constexpr float c_dpt_1 = 0.10;
  constexpr float c_ptmax_1 = 5.0;
  constexpr float c_dpt_2 = 0.20;
  constexpr float c_ptmax_2 = 10.0;
  constexpr float c_dpt_3 = 0.25;
  constexpr float c_dzmax_brl = 0.005;
  constexpr float c_drmax_brl = 0.010;
  constexpr float c_ptmin_hpt = 2.0;
  constexpr float c_dzmax_hpt = 0.010;
  constexpr float c_drmax_hpt = 0.010;
  constexpr float c_dzmax_els = 0.015;
  constexpr float c_drmax_els = 0.015;

  // Threading
  extern int    numThreadsFinder;
  extern int    numThreadsSimulation;

  // For GPU computations
  extern int    numThreadsEvents;
  extern int    numThreadsReorg;

  extern int    finderReportBestOutOfN;

  extern int    numSeedsPerTask;

  // number of layer1 hits for finding seeds per task
  extern int    numHitsPerTask;

  // seed options
  extern seedOpts  seedInput;
  extern cleanOpts seedCleaning;
  
  extern bool   useCMSGeom;
  extern bool   readCmsswTracks;

  extern bool   dumpForPlots;
  extern bool   silent;

  extern bool   kludgeCmsHitErrors;
  extern bool   backwardFit;

  void RecalculateDependentConstants();
  
  CUDA_CALLABLE
  inline float BfieldFromZR(const float z, const float r)
  {
    return (Config::mag_b0*z*z + Config::mag_b1*z + Config::mag_c1)*(Config::mag_a*r*r + 1.f);
  }

#ifdef USE_MATRIPLEX

  #ifndef MPT_SIZE
    #if defined(__MIC__) || defined(__AVX512F__)
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

};

#endif 
