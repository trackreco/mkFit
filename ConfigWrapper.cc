#include "ConfigWrapper.h"
#include "Config.h"
#include "MaterialEffects.h"
#include "TrackerInfo.h"

namespace mkfit {
  namespace ConfigWrapper {
    void initializeForCMSSW(bool silent) {
      Config::seedInput = cmsswSeeds;
      Config::silent = silent;

      // to do backward fit to the first layer, not point of closest approach
      Config::includePCA = false;

      fillZRgridME();
    }

    void setNTotalLayers(int nTotalLayers) {
      Config::nTotalLayers = nTotalLayers;
    }
  }
}
