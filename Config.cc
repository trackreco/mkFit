// Include of Config.h not yet required

namespace Config
{

#ifdef SUPERDEBUG // to avoid include of Config.h
  int nTracks = 1;
  int nEvents = 100000;
#else
  int nTracks = 20000;
  int nEvents = 10;
#endif

  // Multi threading and Clone engine configuration
  int   numThreadsFinder = 1;

#ifdef __MIC__
  int   numThreadsSimulation = 60;
#else
  int   numThreadsSimulation = 12;
#endif

  bool  clonerUseSingleThread  = false;
  int   finderReportBestOutOfN = 1;

  bool  useCMSGeom = false;
}
