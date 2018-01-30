#include "Config.h"

#include "TrackerInfo.h"

namespace Config
{
  TrackerInfo TrkInfo;

  int nTracks = 10000;
  int nEvents = 20;

  int nTotalLayers = -1;

  std::string geomPlugin = "CylCowWLids";

  int maxCandsPerSeed  = 6; // cmssw tests: 6 (GC had 3) \_ set from geom plugin
  int maxHolesPerCand  = 2; // cmssw tests: 12           /

  int maxCandsPerEtaBin; // Recalculated after config is read ... should be removed.

  // Multi threading and Clone engine configuration
  int   numThreadsFinder = 1;
  
  // GPU computations
  int   numThreadsEvents = 1;
  int   numThreadsReorg = 1;

#if defined(__MIC__) || defined(__AVX512F__)
  int   numThreadsSimulation = 60;
#else
  int   numThreadsSimulation = 12;
#endif

  int   finderReportBestOutOfN = 1;

  int   nlayers_per_seed = 3; // can be overriden from Geom plugin; a very confusing variable :)
  int   numSeedsPerTask = 32;
  
  // number of hits per task for finding seeds
  int   numHitsPerTask = 32;

  // material effects
  float RlgridME[Config::nBinsZME][Config::nBinsRME];
  float XigridME[Config::nBinsZME][Config::nBinsRME];

  float chi2Cut = 15.;

  seedOpts  seedInput    = simSeeds;
  cleanOpts seedCleaning = noCleaning; 

  bool             finding_requires_propagation_to_hit_pos;
  PropagationFlags finding_inter_layer_pflags;
  PropagationFlags finding_intra_layer_pflags;
  PropagationFlags backward_fit_pflags;
  PropagationFlags forward_fit_pflags;
  PropagationFlags seed_fit_pflags;
  PropagationFlags pca_prop_pflags;

  bool  usePhiQArrays = true;

  bool  useCMSGeom = false;
  bool  readCmsswTracks = false;

  bool  dumpForPlots = false;
  bool  silent       = false;

  bool  cf_seeding = false;
  bool  cf_fitting = false;

  bool  quality_val = false;
  bool  root_val    = false;
  bool  cmssw_val   = false;
  bool  fit_val     = false;
  bool  readSimTrackStates = false;
  bool  inclusiveShorts = false;
  bool  applyCMSSWHitMatch = false;
  matchOpts cmsswMatching = hitBased;

  bool  kludgeCmsHitErrors = false;
  bool  backwardFit = false;

  void RecalculateDependentConstants()
  {
    maxCandsPerEtaBin = std::max(100, maxCandsPerSeed * (nTracks+100) / nEtaPart);
  }
}
