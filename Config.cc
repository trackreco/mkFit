// Include of Config.h not yet required

namespace Config
{
  // Multi threading and Clone engine configuration
  int   numThreadsFinder = 1;

#ifdef __MIC__
  int   numThreadsSimulation = 60;
#else
  int   numThreadsSimulation = 12;
#endif

  bool  clonerUseSingleThread  = false;
  int   finderReportBestOutOfN = 1;
}
