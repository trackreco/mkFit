//-------------------
// CMS 2017 geometry
//-------------------

#include "Config.h"
#include "TrackerInfo.h"

using namespace mkfit;

namespace
{
#include "CMS-2017.acc"

  void Create_CMS_2017(TrackerInfo& ti, bool verbose)
  {
    Config::nTotalLayers     = 18 + 2 * 27;

    Config::useCMSGeom       = true;
    Config::nlayers_per_seed = 4;
    Config::maxCandsPerSeed  = 6;  // GC said 3 is enough ???
    Config::maxHolesPerCand  = 3; // should be reduced
    Config::chi2Cut          = 15.0;

    Config::finding_requires_propagation_to_hit_pos = true;
    Config::finding_inter_layer_pflags = PropagationFlags(PF_apply_material);
    Config::finding_intra_layer_pflags = PropagationFlags(PF_none);
    Config::backward_fit_pflags        = PropagationFlags(PF_apply_material);
    Config::forward_fit_pflags         = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    Config::seed_fit_pflags            = PropagationFlags(PF_none);
    Config::pca_prop_pflags            = PropagationFlags(PF_none);
    
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
