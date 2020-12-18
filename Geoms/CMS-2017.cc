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

  void IterationConfig::import_seeds(Event *ev, TrackerInfo& trk_info, unsigned int it=0){
    
    // Seeds are placed into eta regions and sorted on eta. Counts for each eta region are
    // stored into Event::seedEtaSeparators_.
    
    //bool debug = true;
    
    MkSeedPacket mkSeedPacket_(ev);
    
    //////TrackerInfo &trk_info = Config::TrkInfo;
    TrackVec    &seeds    = mkSeedPacket_.m_inseeds_; //
    const int    size     = seeds.size();

    for (int i = 0; i < 5; ++i)
      {
	mkSeedPacket_.m_seedEtaSeparators_[i] = 0; // To mkbuilder
	mkSeedPacket_.m_seedMinLastLayer_ [i] = 9999;
	mkSeedPacket_.m_seedMaxLastLayer_ [i] = 0;
      }
    
    std::vector<float> etas(size);
    for (int i = 0; i < size; ++i)
      {
	const Track &S         = seeds[i];
	const bool   z_dir_pos = S.pz() > 0;
	
	HitOnTrack hot = S.getLastHitOnTrack();
	float      eta = mkSeedPacket_.m_layerHits_[hot.layer][hot.index].eta();
	// float   eta = S.momEta();
	
	// Region to be defined by propagation / intersection tests
	TrackerInfo::EtaRegion reg;
	
	// Hardcoded for cms ... needs some lists of layers (hit/miss) for brl / ecp tests.
	// MM: Check lambda functions/std::function 
	const LayerInfo &outer_brl = trk_info.outer_barrel_layer();
	
	const LayerInfo &tib1 = trk_info.m_layers[ 4];
	const LayerInfo &tob1 = trk_info.m_layers[10];
	
	const LayerInfo &tecp1 = trk_info.m_layers[27];
	const LayerInfo &tecn1 = trk_info.m_layers[54];
	
	const LayerInfo &tec_first = z_dir_pos ? tecp1 : tecn1;
	
	// If a track hits outer barrel ... it is in the barrel (for central, "outgoing" tracks).
	// This is also true for cyl-cow.
	// Better check is: hits outer TIB, misses inner TEC (but is +-z dependant).
	// XXXX Calculate z ... then check is inside or less that first EC z.
	// There are a lot of tracks that go through that crack.
	
	// XXXX trying a fix for low pT tracks that are in barrel after half circle
	float maxR = S.maxReachRadius();
	float z_at_maxr;
	
	bool  can_reach_outer_brl = S.canReachRadius(outer_brl.m_rout);
	float z_at_outer_brl;
	bool  misses_first_tec;
	if (can_reach_outer_brl)
	  {
	    z_at_outer_brl = S.zAtR(outer_brl.m_rout);
	    if (z_dir_pos)
	      misses_first_tec = z_at_outer_brl < tec_first.m_zmin;
	    else
	      misses_first_tec = z_at_outer_brl > tec_first.m_zmax;
	  }
	else
	  {
	    z_at_maxr = S.zAtR(maxR);
	    if (z_dir_pos)
	      misses_first_tec = z_at_maxr < tec_first.m_zmin;
	    else
	      misses_first_tec = z_at_maxr > tec_first.m_zmax;
	  }
	
	if (/*can_reach_outer_brl &&*/ misses_first_tec)
	  // outer_brl.is_within_z_limits(S.zAtR(outer_brl.r_mean())))
	  {
	    reg = TrackerInfo::Reg_Barrel;
	  }
	else
	  {
	    // This should be a list of layers
	    // CMS, first tib, tob: 4, 10
	    
	    if ((S.canReachRadius(tib1.m_rin) && tib1.is_within_z_limits(S.zAtR(tib1.m_rin))) ||
		(S.canReachRadius(tob1.m_rin) && tob1.is_within_z_limits(S.zAtR(tob1.m_rin))) )
	      {
		// transition region ... we are still hitting barrel layers
		
		reg = z_dir_pos ? TrackerInfo::Reg_Transition_Pos : TrackerInfo::Reg_Transition_Neg;
	      }
	    else
	      {
		// endcap ... no barrel layers will be hit anymore.
		
		reg = z_dir_pos ? TrackerInfo::Reg_Endcap_Pos : TrackerInfo::Reg_Endcap_Neg;
	      }
	  }
	// reg is now defined
	
	++mkSeedPacket_.m_seedEtaSeparators_[reg];
	
	mkSeedPacket_.m_seedMinLastLayer_[reg] = std::min(mkSeedPacket_.m_seedMinLastLayer_[reg], hot.layer);
	mkSeedPacket_.m_seedMaxLastLayer_[reg] = std::max(mkSeedPacket_.m_seedMaxLastLayer_[reg], hot.layer);
	
	etas[i] = 5.0f * (reg - 2) + eta;
	
      }
    
    for (int i = 0; i < 5; ++i)
      {
	if (mkSeedPacket_.m_seedMinLastLayer_[i] == 9999) mkSeedPacket_.m_seedMinLastLayer_[i] = -1;
	if (mkSeedPacket_.m_seedMaxLastLayer_[i] ==    0) mkSeedPacket_.m_seedMaxLastLayer_[i] = -1;
      }
    
    RadixSort rs;
    rs.Sort(&etas[0], size);
    
    TrackVec orig_seeds;
    orig_seeds.swap(seeds);
    seeds.reserve(size);
    for (int i = 0; i < size; ++i)
      {
	seeds.emplace_back( orig_seeds[ rs.GetRanks()[i] ] );
      }
    
    dprintf("IterationConfig::import_seeds finished import of %d seeds (last seeding layer min, max):\n"
	    "  ec- = %d(%d,%d), t- = %d(%d,%d), brl = %d(%d,%d), t+ = %d(%d,%d), ec+ = %d(%d,%d).\n",
	    size,
	    mkSeedPacket_.m_seedEtaSeparators_[0], mkSeedPacket_.m_seedMinLastLayer_[0], mkSeedPacket_.m_seedMaxLastLayer_[0],
	    mkSeedPacket_.m_seedEtaSeparators_[1], mkSeedPacket_.m_seedMinLastLayer_[1], mkSeedPacket_.m_seedMaxLastLayer_[1],
	    mkSeedPacket_.m_seedEtaSeparators_[2], mkSeedPacket_.m_seedMinLastLayer_[2], mkSeedPacket_.m_seedMaxLastLayer_[2],
	    mkSeedPacket_.m_seedEtaSeparators_[3], mkSeedPacket_.m_seedMinLastLayer_[3], mkSeedPacket_.m_seedMaxLastLayer_[3],
	    mkSeedPacket_.m_seedEtaSeparators_[4], mkSeedPacket_.m_seedMinLastLayer_[4], mkSeedPacket_.m_seedMaxLastLayer_[4]);
    
    // Sum region counts up to contain actual separator indices:
    for (int i = TrackerInfo::Reg_Transition_Neg; i < TrackerInfo::Reg_Count; ++i)
      {
	mkSeedPacket_.m_seedEtaSeparators_[i] += mkSeedPacket_.m_seedEtaSeparators_[i - 1];
      }
    
  }
  
  void Create_IterationConfig(Event *ev, TrackerInfo& ti, IterationParams& ip, unsigned int it=0){
    
    thisIter = IterationConfig::IterationConfig(ti, ip, it);
    setIterSteeringParams                      (ti, ip, it);
    setIterRegions                             (ti, ip, it);
    setIterParams                              (    ip, it);
    import_seeds                               (ev, ti, it);
  }

  void Create_CMS_2017(Event *ev, TrackerInfo& ti, IterationParams& ip, unsigned int iter=0, bool verbose)
  {
    
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

    if (iter==0)
      Create_IterationConfig(ev, ti, ip, iter);

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
