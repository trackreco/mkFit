#include "ConfigWrapper.h"
#include "Config.h"
#include "MaterialEffects.h"
#include "TrackerInfo.h"

namespace mkfit {
  namespace ConfigWrapper {
    void initializeForCMSSW(SeedCleaningOpts seedClean, BackwardFit backfit, bool silent) {
      Config::seedInput = cmsswSeeds;
      Config::geomPlugin = "CMS-2017";
      Config::silent = silent;
      Config::cmssw_export = true;

      if(seedClean == SeedCleaningOpts::cleanSeedsN2) {
        Config::seedCleaning = cleanSeedsN2;
      }

      switch(backfit) {
      case BackwardFit::noFit:
        Config::backwardFit = false;
        break;
      case BackwardFit::toFirstLayer:
        Config::backwardFit = true;
        Config::includePCA = false;
        break;
      case BackwardFit::toPCA:
        Config::backwardFit = true;
        Config::includePCA = true;
        break;
      }

      TrackerInfo::ExecTrackerInfoCreatorPlugin(Config::geomPlugin, Config::TrkInfo);

      fillZRgridME();
    }

    void setRemoveDuplicates(bool removeDuplicates) {
      Config::removeDuplicates = removeDuplicates;
    }

    void setNTotalLayers(int nTotalLayers) {
      Config::nTotalLayers = nTotalLayers;
    }
  }
}
