//-------------------
// CMS 2017 geometry
//-------------------

#include "Config.h"
#include "TrackerInfo.h"

namespace
{
#include "CMS-2017.acc"

  void Create_CMS_2017(TrackerInfo& ti, bool verbose)
  {
    Config::nTotalLayers     = 18 + 2 * 27;
    Config::useCMSGeom       = true;
    Config::nlayers_per_seed = 4;
    Config::maxCandsPerSeed  = 6;  // GC said 3 is enough ???
    Config::maxHolesPerCand  = 12; // should be reduced
    Config::chi2Cut          = 30.0;

    ti.set_eta_regions(0.9, 1.7, 2.45, false);
    ti.create_layers(18, 27, 27);

    Create_CMS_2017_AutoGen(ti);

    if (verbose)
    {
      printf("==========================================================================================\n");
    }

    printf("CMS-2017 -- Create_TrackerInfo finished\n");

    if (verbose)
    {
      printf("==========================================================================================\n");
      for (auto &i : ti.m_layers)  i.print_layer();
      printf("==========================================================================================\n");
    }
  }
}


void* TrackerInfoCrator_ptr = (void*) Create_CMS_2017;
