//-------------------
// CMS 2017 geometry
//-------------------

#include "Config.h"
#include "TrackerInfo.h"
#include "mkFit/SteeringParams.h"

using namespace mkfit;

namespace
{
#include "CMS-2017.acc"

  void setIterSteeringParams(TrackerInfo& ti, IterationParams& ip, unsigned int it=0)
  {
    if(it==0)      
      {	
	{ SteeringParams &sp = m_steering_params[TrackerInfo::Reg_Endcap_Neg];
	  sp.reserve_plan(3 + 3 + 6 + 18);
	  sp.fill_plan(0, 1, false, true); // bk-fit only
	  sp.append_plan( 2, true);        // pick-up only
	  sp.append_plan(45, false);
	  sp.append_plan(46, false);
	  sp.append_plan(47, false);
	  sp.fill_plan(48, 53); // TID,  6 layers
	  sp.fill_plan(54, 71); // TEC, 18 layers
	  sp.finalize_plan();
	}
	
	{ SteeringParams &sp = m_steering_params[TrackerInfo::Reg_Transition_Neg];
	  sp.reserve_plan(3 + 4 + 6 + 6 + 8 + 18);
	  sp.fill_plan(0, 1, false, true); // bk-fit only
	  sp.append_plan( 2, true);
	  sp.append_plan( 3, false);
	  sp.append_plan(45, false);
	  sp.append_plan(46, false);
	  sp.append_plan(47, false);
	  sp.fill_plan( 4,  9); // TIB,  6 layers
	  sp.fill_plan(48, 53); // TID,  6 layers
	  sp.fill_plan(10, 17); // TOB,  8 layers
	  sp.fill_plan(54, 71); // TEC, 18 layers
	  sp.finalize_plan();
	}
	
	{ SteeringParams &sp = m_steering_params[TrackerInfo::Reg_Barrel];
	  sp.reserve_plan(3 + 1 + 6 + 8);
	  sp.fill_plan(0, 1, false, true); // bk-fit only
	  sp.append_plan( 2, true);        // pickup-only
	  sp.append_plan( 3, false);
	  sp.fill_plan( 4,  9); // TIB, 6 layers
	  sp.fill_plan(10, 17); // TOB, 8 layers
	  sp.finalize_plan();
	}
	
	{ SteeringParams &sp = m_steering_params[TrackerInfo::Reg_Transition_Pos];
	  sp.reserve_plan(3 + 4 + 6 + 6 + 8 + 18);
	  sp.fill_plan(0, 1, false, true); // bk-fit only
	  sp.append_plan( 2, true);        // pickup-only
	  sp.append_plan( 3, false);
	  sp.append_plan(18, false);
	  sp.append_plan(19, false);
	  sp.append_plan(20, false);
	  sp.fill_plan( 4,  9); // TIB,  6 layers
	  sp.fill_plan(21, 26); // TID,  6 layers
	  sp.fill_plan(10, 17); // TOB,  8 layers
	  sp.fill_plan(27, 44); // TEC, 18 layers
	  sp.finalize_plan();
	}
	
	{ SteeringParams &sp = m_steering_params[TrackerInfo::Reg_Endcap_Pos];
	  sp.reserve_plan(3 + 3 + 6 + 18);
	  sp.fill_plan(0, 1, false, true); // bk-fit only
	  sp.append_plan( 2, true);        // pickup-only
	  sp.append_plan(18, false);
	  sp.append_plan(19, false);
	  sp.append_plan(20, false);
	  sp.fill_plan(21, 26); // TID,  6 layers
	  sp.fill_plan(27, 44); // TEC, 18 layers
	  sp.finalize_plan();
	}	

      }
  }

  void setIterRegions(TrackerInfo& ti, IterationParams& ip, unsigned int it=0)
  {  
    if(it==0)
      {
      m_regions.resize(5);
      m_regions[0] = TrackerInfo::Reg_Transition_Pos;
      m_regions[1] = TrackerInfo::Reg_Transition_Neg;
      m_regions[2] = TrackerInfo::Reg_Endcap_Pos;
      m_regions[3] = TrackerInfo::Reg_Endcap_Neg;
      m_regions[4] = TrackerInfo::Reg_Barrel;
      }
  }

  void setIterationParams(IterationParams& ip, unsigned int it=0)
  {
    if(it==0)
      {
	nlayers_per_seed = 4;
	maxCandsPerSeed  = 5;
	maxHolesPerCand  = 4;
	maxConsecHoles   = 1;
	chi2Cut          = 30;
      }
  }

  void Create_IterationConfig(TrackerInfo& ti, IterationParams& ip, unsigned int it){
    
    thisIter = IterationConfig::IterationConfig(ti, ip, it);
    setIterSteeringParams                      (ti, ip, it);
    setIterRegions                             (ti, ip, it);
    setIterParams                              (    ip, it);

  }

  void Create_CMS_2017(TrackerInfo& ti, IterationParams& ip, unsigned int iter=0, bool verbose)
  {
    
    if (iter==0)
      Create_IterationConfig(ti, ip, iter);

    Config::nTotalLayers     = 18 + 2 * 27;

    Config::useCMSGeom       = true;

// MIMI -- Pandora's Box -- to be commented out, moved to iteration config.
    Config::maxCandsPerSeed  = 5;
    Config::maxHolesPerCand  = 4;
    Config::maxConsecHoles   = 1;
    Config::chi2Cut          = 30;
    Config::chi2CutOverlap   = 3.5;
    Config::pTCutOverlap     = 1;

// MIMI -- This was from Mario:
//The following are commented out, since now set by IterationConfig (once per iteration)
//    Config::nlayers_per_seed = 4;
//    Config::maxCandsPerSeed  = 5;
//    Config::maxHolesPerCand  = 4;
//    Config::maxConsecHoles   = 1;
//    Config::chi2Cut          = 30;

    Config::finding_requires_propagation_to_hit_pos = true;
    Config::finding_inter_layer_pflags = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    Config::finding_intra_layer_pflags = PropagationFlags(PF_none);
    Config::backward_fit_pflags        = PropagationFlags(PF_use_param_b_field | PF_apply_material);
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
