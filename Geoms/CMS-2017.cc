//-------------------
// CMS 2017 geometry
//-------------------

#include "../Config.h"
#include "../TrackerInfo.h"

namespace
{
#include "CMS-2017.acc"

  void Create_CMS_2017(TrackerInfo& ti, bool verbose)
  {
    Config::useCMSGeom = true;
    Config::chi2Cut    = 30.0;

    ti.set_eta_regions(0.9, 1.7, 2.45);
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
