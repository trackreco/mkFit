#include "Config.h"

#include "TrackerInfo.h"

namespace Config
{
  TrackerInfo TrkInfo;

  int nTracks = 10000;
  int nEvents = 20;

  std::string geomPlugin = "CylCowWLids";

  // Dependent constants, assigned after processing of commandline options
  int maxHitsPerBunch;
  int maxCandsPerEtaBin;

  // Multi threading and Clone engine configuration
  int   numThreadsFinder = 1;
  
  // GPU computations
  int   numThreadsEvents = 1;
  int   numThreadsReorg = 1;

#ifdef __MIC__
  int   numThreadsSimulation = 60;
#else
  int   numThreadsSimulation = 12;
#endif

  int   finderReportBestOutOfN = 1;

  int   nlayers_per_seed = 3; // default is 3 for barrel seeding --> will need a new variable once we move to endcap seeding
  int   numSeedsPerTask = 32;
  
  // number of hits per task for finding seeds
  int   numHitsPerTask = 32;

  // material effects
  float RlgridME[Config::nBinsZME][Config::nBinsRME];
  float XigridME[Config::nBinsZME][Config::nBinsRME];

  float chi2Cut = 15.;

  bool  useCMSGeom = false;
  bool  readCmsswSeeds = false;

  bool  findSeeds  = false;
  bool  endcapTest = false;

  bool  silent     = false;

  bool  cf_seeding = false;
  bool  cf_fitting = false;

  bool  root_val = false;
  bool  fit_val  = false;
  bool  inclusiveShorts = false;

  void RecalculateDependentConstants()
  {
    maxCandsPerEtaBin = std::max(100, maxCandsPerSeed * (nTracks+100) / nEtaPart);
    maxHitsPerBunch   = std::max(100, nTracks * 12 / 10 / nEtaPart) + maxHitsConsidered;
  }
}
