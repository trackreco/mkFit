#include "Config.h"

#include "TrackerInfo.h"

namespace mkfit {

namespace Config
{
  TrackerInfo TrkInfo;

  int nTracks = 10000;
  int nEvents = 20;
  bool loopOverFile = false;

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

#ifdef CONFIG_PhiQArrays
  bool  usePhiQArrays = true;
#endif

  bool  useCMSGeom = false;
  bool  readCmsswTracks = false;

  bool  dumpForPlots = false;
  bool  silent       = false;

  bool  cf_seeding = false;
  bool  cf_fitting = false;

  bool  quality_val = false;
  bool  sim_val_for_cmssw = false;
  bool  sim_val     = false;
  bool  cmssw_val   = false;
  bool  cmssw_export = false;
  bool  fit_val     = false;
  bool  readSimTrackStates = false;
  bool  inclusiveShorts = false;
  bool  keepHitInfo = false;
  bool  tryToSaveSimInfo = false;
  matchOpts cmsswMatchingFW = hitBased;
  matchOpts cmsswMatchingBK = trkParamBased;

  bool  removeDuplicates = false;
  bool  useHitsForDuplicates = true;
  float maxdPhi = 0.5;
  float maxdPt  = 0.5;
  float maxdEta = 0.05;
  float minFracHitsShared = 0.75;
  float maxdRSquared = 0.000001; //corresponds to maxdR of 0.001

  bool mtvLikeValidation = false;
  bool mtvRequireSeeds = false;
  int  cmsSelMinLayers = 12;
  int  nMinFoundHits = 10;

  bool  kludgeCmsHitErrors = false;
  bool  backwardFit = false;
  bool  includePCA = false;

  void RecalculateDependentConstants()
  {
    maxCandsPerEtaBin = std::max(100, maxCandsPerSeed * (nTracks+100) / nEtaPart);
  }
}

} // end namespace mkfit
