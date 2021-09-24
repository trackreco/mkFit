//-------------------
// CMS 2017 geometry
//-------------------

#include "Config.h"
#include "Debug.h"
#include "TrackerInfo.h"
#include "mkFit/IterationConfig.h"
#include "mkFit/HitStructures.h"

#include "CMS-2017-HitSelectionWindows.h"

#include <functional>

using namespace mkfit;

namespace
{
#include "CMS-2017.acc"

  void SetupCoreSteeringParams_Iter0(IterationConfig& ic)
  {
    ic.m_region_order[0] = TrackerInfo::Reg_Transition_Pos;
    ic.m_region_order[1] = TrackerInfo::Reg_Transition_Neg;
    ic.m_region_order[2] = TrackerInfo::Reg_Endcap_Pos;
    ic.m_region_order[3] = TrackerInfo::Reg_Endcap_Neg;
    ic.m_region_order[4] = TrackerInfo::Reg_Barrel;

    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Endcap_Neg];
      sp.reserve_plan(3 + 3 + 6 + 18); // BPix + FPix- + TID- + TEC-; BPix4 is out of acceptance
      sp.fill_plan( 0,  2);
      sp.fill_plan(45, 47);
      sp.fill_plan(48, 53); // TID,  6 disks (3 mono + 3 stereo)
      sp.fill_plan(54, 71); // TEC, 18 disks (3 mono + 3 stereo)
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Transition_Neg];
      sp.reserve_plan(3 + 4 + 6 + 6 + 8 + 18); // BPix + FPix- + TIB + TID- + TOB + TEC-
      sp.fill_plan (0,  3);
      sp.fill_plan(45, 47);
      sp.fill_plan( 4,  9); // TIB,  6 layers (4 mono + 2 stereo)
      sp.fill_plan(48, 53); // TID,  6 disks  (3 mono + 3 stereo)
      sp.fill_plan(10, 17); // TOB,  8 layers (6 mono + 2 stereo)
      sp.fill_plan(54, 71); // TEC, 18 disks  (9 mono + 9 stereo)
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Barrel];
      sp.reserve_plan(4 + 6 + 8); // BPix + TIB + TOB
      sp.fill_plan( 0,  3);
      sp.fill_plan( 4,  9); // TIB, 6 layers (4 mono + 2 stereo)  [ 4,  9]
      sp.fill_plan(10, 17); // TOB, 8 layers (6 mono + 2 stereo)  [10, 17]
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Transition_Pos];
      sp.reserve_plan(3 + 4 + 6 + 6 + 8 + 18); // BPix + FPix+ + TIB + TID+ + TOB + TEC+
      sp.fill_plan( 0,  3);
      sp.fill_plan(18, 20);
      sp.fill_plan( 4,  9); // TIB,  6 layers (4 mono + 2 stereo)  [ 7, 12]
      sp.fill_plan(21, 26); // TID,  6 disks  (3 mono + 3 stereo)  [13, 18]
      sp.fill_plan(10, 17); // TOB,  8 layers (6 mono + 2 stereo)  [19, 26]
      sp.fill_plan(27, 44); // TEC, 18 disks  (9 mono + 9 stereo)  [27, 44]
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Endcap_Pos];
      sp.reserve_plan(3 + 3 + 6 + 18); // BPix + FPix+ + TID+ + TEC+; BPix4 is out of acceptance
      sp.fill_plan( 0,  2);
      sp.fill_plan(18, 20);
      sp.fill_plan(21, 26); // TID,  6 disks  (3 mono + 3 stereo)  [ 6, 11]
      sp.fill_plan(27, 44); // TEC, 18 disks  (9 mono + 9 stereo)  [12, 29]
      sp.set_iterator_limits(2, 0);
    }
  }

  void SetupBackwardSearch_Iter0(IterationConfig& ic)
  {
    ic.m_backward_params = ic.m_params;
    ic.m_backward_search = true;
    ic.m_backward_drop_seed_hits = true;
    auto &spv = ic.m_steering_params;
    spv[TrackerInfo::Reg_Endcap_Neg]    .set_iterator_limits(2, 3, 5);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(2, 3, 7);
    spv[TrackerInfo::Reg_Barrel]        .set_iterator_limits(2, 3, 4);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(2, 3, 7);
    spv[TrackerInfo::Reg_Endcap_Pos]    .set_iterator_limits(2, 3, 5);
  }

  void SetupBackwardSearch_Iter5(IterationConfig& ic)
  {
    ic.m_backward_params = ic.m_params;
    ic.m_backward_search = true;
    ic.m_backward_drop_seed_hits = true;
    auto &spv = ic.m_steering_params;
    spv[TrackerInfo::Reg_Endcap_Neg]    .set_iterator_limits(1, 2, 7);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(1, 2, 9);
    spv[TrackerInfo::Reg_Barrel]        .set_iterator_limits(1, 2, 6);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(1, 2, 9);
    spv[TrackerInfo::Reg_Endcap_Pos]    .set_iterator_limits(1, 2, 7);
  }

  void SetupBackwardSearch_Iter7(IterationConfig& ic)
  {
    ic.m_backward_search = true;
    ic.m_backward_params = ic.m_params;
    ic.m_backward_params.maxHolesPerCand = 2;
    ic.m_backward_params.maxConsecHoles  = 2;
    // Remove pixel layers from FwdSearch, add them to BkwSearch
    auto &spv = ic.m_steering_params;
    spv[TrackerInfo::Reg_Endcap_Neg]    .set_iterator_limits(8, 6, 19);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(9, 7, 34);
    spv[TrackerInfo::Reg_Barrel]        .set_iterator_limits(6, 4, 8);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(9, 7, 34);
    spv[TrackerInfo::Reg_Endcap_Pos]    .set_iterator_limits(8, 6, 19);
  }

  void SetupBackwardSearch_Iter8(IterationConfig& ic)
  {
    ic.m_backward_search = true;
    ic.m_backward_params = ic.m_params;
    ic.m_backward_params.maxHolesPerCand = 2;
    ic.m_backward_params.maxConsecHoles  = 2;
    // Remove pixel/tib/tid layers from FwdSearch, add them to BkwSearch/
    auto &spv = ic.m_steering_params;
    spv[TrackerInfo::Reg_Endcap_Neg]    .set_iterator_limits(12, 12, 24);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(27, 19, 39);
    spv[TrackerInfo::Reg_Barrel]        .set_iterator_limits(12, 10, 14);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(27, 19, 39);
    spv[TrackerInfo::Reg_Endcap_Pos]    .set_iterator_limits(12, 12, 24);
  }

  void SetupIterationParams(IterationParams& ip, unsigned int it=0)
  {
    if (it == 0)
    {
      ip.nlayers_per_seed = 4;
      ip.maxCandsPerSeed  = 5;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 10.0;
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    }
    else if (it == 1) // for triplet steps, nlayers_per_seed=3
    {   
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed  = 5;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 10.0; 
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    }
    else if (it == 2)
    {   
      ip.nlayers_per_seed = 4;
      ip.maxCandsPerSeed  = 5;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 10.0; 
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    } 
    else if (it == 3) // for triplet steps, nlayers_per_seed=3
    {   
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed  = 5;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 10.0; 
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    } 
    else if (it == 4)
    {   
      ip.nlayers_per_seed = 4;
      ip.maxCandsPerSeed  = 5;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 10.0; 
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    } 
    else if (it == 5) // for triplet steps, nlayers_per_seed=3
    {   
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed  = 5;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 10.0; 
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    } 
    else if (it == 6) // for triplet steps, nlayers_per_seed=3; for mixeTripletSetp, also maxCandsPerSeed=2
    {
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed  = 2;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 10.0;
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    } 
    else if (it == 7) // for PixelLess step, maxCandsPerSeed=2 and maxHolesPerCand=maxConsecHoles=0
    {   
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed  = 2;
      ip.maxHolesPerCand  = 0;
      ip.maxConsecHoles   = 1;
      ip.chi2Cut_min      = 10.0; 
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    }   
    else if (it == 8) // for TobTec step, maxCandsPerSeed=2 and maxHolesPerCand=maxConsecHoles=0
    {   
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed  = 2;
      ip.maxHolesPerCand  = 0;
      ip.maxConsecHoles   = 1;
      ip.chi2Cut_min      = 10.0; 
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    }
    else if (it == 9) // addign also pixel pair step - algo -> 6
    {
      ip.nlayers_per_seed = 2;
      ip.maxCandsPerSeed  = 3;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 10.0;
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    }
  }

  void fill_hit_selection_windows_params(IterationConfig &ic)
  {
    HitSelectionWindows hsw;
    for (int l = 0; l < (int)ic.m_layer_configs.size(); ++l)
    {
      // dphi cut
      ic.m_layer_configs[l].c_dp_0 = hsw.m_dp_params[ic.m_iteration_index][l][0];
      ic.m_layer_configs[l].c_dp_1 = hsw.m_dp_params[ic.m_iteration_index][l][1];
      ic.m_layer_configs[l].c_dp_2 = hsw.m_dp_params[ic.m_iteration_index][l][2];
      // dq cut
      ic.m_layer_configs[l].c_dq_0 = hsw.m_dq_params[ic.m_iteration_index][l][0];
      ic.m_layer_configs[l].c_dq_1 = hsw.m_dq_params[ic.m_iteration_index][l][1];
      ic.m_layer_configs[l].c_dq_2 = hsw.m_dq_params[ic.m_iteration_index][l][2];
      // chi2 cut (for future optimization)
      ic.m_layer_configs[l].c_c2_0 = hsw.m_c2_params[ic.m_iteration_index][l][0];
      ic.m_layer_configs[l].c_c2_1 = hsw.m_c2_params[ic.m_iteration_index][l][1];
      ic.m_layer_configs[l].c_c2_2 = hsw.m_c2_params[ic.m_iteration_index][l][2];
    }
  }

  std::function<IterationConfig::partition_seeds_foo> PartitionSeeds0 =
  [](const TrackerInfo &trk_info, const TrackVec &in_seeds, const EventOfHits &eoh,
     IterationSeedPartition &part)
  {
    // Seeds are placed into eta regions and sorted on region + eta.

    const int size = in_seeds.size();

    for (int i = 0; i < size; ++i)
    {
      const Track &S = in_seeds[i];

      const bool z_dir_pos = S.pz() > 0;

      HitOnTrack hot = S.getLastHitOnTrack();
      // MIMI ACHTUNG -- here we assume seed hits have already been remapped.
      // This was true at that time :)
      float eta = eoh[hot.layer].GetHit(hot.index).eta();
      // float  eta = S.momEta();

      // Region to be defined by propagation / intersection tests
      TrackerInfo::EtaRegion reg;

      // Hardcoded for cms ... needs some lists of layers (hit/miss) for brl / ecp tests.
      // MM: Check lambda functions/std::function
      const LayerInfo &outer_brl = trk_info.outer_barrel_layer();

      const LayerInfo &tib1 = trk_info.m_layers[4];
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
            (S.canReachRadius(tob1.m_rin) && tob1.is_within_z_limits(S.zAtR(tob1.m_rin))))
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

      part.m_region[i]     = reg;
      part.m_sort_score[i] = 5.0f * (reg - 2) + eta;
    }
  };


  void Create_CMS_2017(TrackerInfo& ti, IterationsInfo& ii, bool verbose)
  {
    Config::nTotalLayers     = 18 + 2 * 27;

    Config::useCMSGeom       = true;

    Config::finding_requires_propagation_to_hit_pos = true;
    Config::finding_inter_layer_pflags = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    Config::finding_intra_layer_pflags = PropagationFlags(PF_none);
    Config::backward_fit_pflags        = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    Config::forward_fit_pflags         = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    Config::seed_fit_pflags            = PropagationFlags(PF_none);
    Config::pca_prop_pflags            = PropagationFlags(PF_none);

    ti.set_eta_regions(0.9, 1.7, 2.45, false);
    ti.create_layers(18, 27, 27);

    ii.resize(10);
    ii[0].set_iteration_index_and_track_algorithm(0, (int) TrackBase::TrackAlgorithm::initialStep);
    ii[0].set_num_regions_layers(5, 72);

    // Fills TrackerInfo/LayerInfo and default windows of ii[0].m_layer_configs
    Create_CMS_2017_AutoGen(ti, ii);
    ii[0].m_partition_seeds = PartitionSeeds0;

    SetupCoreSteeringParams_Iter0(ii[0]);

    // At this point copy out layer/steering stuff for reuse in later iterations.
    IterationConfig def_itconf;
    def_itconf.CloneLayerSteerCore(ii[0]);

    SetupIterationParams(ii[0].m_params, 0);
    ii[0].set_dupclean_flag();
    ii[0].set_dupl_params(0.5, 0.002,0.004,0.008);
    fill_hit_selection_windows_params(ii[0]);
    // Backward-search with seed region rebuilding
    SetupBackwardSearch_Iter0(ii[0]);

    ii[1].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[1].m_params, 1);
    ii[1].set_iteration_index_and_track_algorithm(1, (int) TrackBase::TrackAlgorithm::highPtTripletStep);
    ii[1].set_seed_cleaning_params(2.0, 0.018, 0.018, 0.018, 0.018, 0.018, 0.05, 0.018, 0.05); 
    ii[1].set_dupclean_flag();
    ii[1].set_dupl_params(0.5, 0.03,0.05,0.08);
    fill_hit_selection_windows_params(ii[1]);
    ii[1].m_backward_params = ii[1].m_params;

    ii[2].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[2].m_params, 2);
    ii[2].set_iteration_index_and_track_algorithm(2, (int) TrackBase::TrackAlgorithm::lowPtQuadStep);
    ii[2].set_seed_cleaning_params(0.5, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05);
    ii[2].set_dupclean_flag();
    ii[2].set_dupl_params(0.5, 0.01,0.03,0.05);
    fill_hit_selection_windows_params(ii[2]);
    ii[2].m_backward_params = ii[2].m_params;
     
    ii[3].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[3].m_params, 3);
    ii[3].set_iteration_index_and_track_algorithm(3, (int) TrackBase::TrackAlgorithm::lowPtTripletStep);
    ii[3].set_seed_cleaning_params(0.5, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05);
    ii[3].set_dupclean_flag();
    ii[3].set_dupl_params(0.33, 0.018,0.05,0.018);
    fill_hit_selection_windows_params(ii[3]);
    ii[3].m_backward_params = ii[3].m_params;
    
    ii[4].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[4].m_params, 4);
    ii[4].set_iteration_index_and_track_algorithm(4, (int) TrackBase::TrackAlgorithm::detachedQuadStep);
    ii[4].set_seed_cleaning_params(2.0, 0.018, 0.018, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05);
    ii[4].set_dupclean_flag();
    ii[4].set_dupl_params(0.25, 0.018,0.05,0.05);
    fill_hit_selection_windows_params(ii[4]);
    ii[4].m_backward_params = ii[4].m_params;
    
    ii[5].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[5].m_params, 5);
    ii[5].set_iteration_index_and_track_algorithm(5, (int) TrackBase::TrackAlgorithm::detachedTripletStep);
    ii[5].set_seed_cleaning_params(2.0, 0.018, 0.018, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05);
    ii[5].set_dupclean_flag();
    ii[5].set_dupl_params(0.25, 0.05,0.05,0.05);
    fill_hit_selection_windows_params(ii[5]);
    // Backward-search with seed region rebuilding
    SetupBackwardSearch_Iter5(ii[5]);

    ii[6].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[6].m_params, 6);
    ii[6].set_iteration_index_and_track_algorithm(6, (int) TrackBase::TrackAlgorithm::mixedTripletStep);
    ii[6].set_seed_cleaning_params(2.0, 0.05, 0.05, 0.135, 0.135, 0.05, 0.05, 0.135, 0.135);
    ii[6].set_dupclean_flag();
    ii[6].set_dupl_params(0.2, 0.05,0.05,0.05); 
    fill_hit_selection_windows_params(ii[6]);
    ii[6].m_backward_params = ii[6].m_params;

    ii[7].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[7].m_params, 7);
    ii[7].set_iteration_index_and_track_algorithm(7, (int) TrackBase::TrackAlgorithm::pixelLessStep);
    ii[7].set_seed_cleaning_params(2.0, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135);
    ii[7].set_qf_flags();
    ii[7].set_qf_params(4,0.19);
    fill_hit_selection_windows_params(ii[7]);
    SetupBackwardSearch_Iter7(ii[7]);

    ii[8].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[8].m_params, 8);
    ii[8].set_iteration_index_and_track_algorithm(8, (int) TrackBase::TrackAlgorithm::tobTecStep);
    ii[8].set_seed_cleaning_params(2.0, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135);    
    ii[8].set_qf_flags();
    ii[8].set_qf_params(4,0.25);
    fill_hit_selection_windows_params(ii[8]);
    SetupBackwardSearch_Iter8(ii[8]);

    ii[9].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[9].m_params, 9);
    ii[9].set_iteration_index_and_track_algorithm(9, (int) TrackBase::TrackAlgorithm::pixelPairStep);
    ii[9].set_seed_cleaning_params(2.0, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135);
    ii[9].set_dupclean_flag();
    ii[9].set_dupl_params(0.5, 0.03,0.05,0.05);
    fill_hit_selection_windows_params(ii[9]);
    ii[9].m_backward_params = ii[9].m_params;

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
