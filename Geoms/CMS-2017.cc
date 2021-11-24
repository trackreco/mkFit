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
      sp.fill_plan( 0,  3); //                                    [ 0,  3]
      sp.fill_plan( 4,  9); // TIB, 6 layers (4 mono + 2 stereo)  [ 4,  9]
      sp.fill_plan(10, 17); // TOB, 8 layers (6 mono + 2 stereo)  [10, 17]
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Transition_Pos];
      sp.reserve_plan(3 + 4 + 6 + 6 + 8 + 18); // BPix + FPix+ + TIB + TID+ + TOB + TEC+
      sp.fill_plan( 0,  3); //                                     [ 0,  3]
      sp.fill_plan(18, 20); //                                     [ 4,  6]
      sp.fill_plan( 4,  9); // TIB,  6 layers (4 mono + 2 stereo)  [ 7, 12]
      sp.fill_plan(21, 26); // TID,  6 disks  (3 mono + 3 stereo)  [13, 18]
      sp.fill_plan(10, 17); // TOB,  8 layers (6 mono + 2 stereo)  [19, 26]
      sp.fill_plan(27, 44); // TEC, 18 disks  (9 mono + 9 stereo)  [27, 44]
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Endcap_Pos];
      sp.reserve_plan(3 + 3 + 6 + 18); // BPix + FPix+ + TID+ + TEC+; BPix4 is out of acceptance
      sp.fill_plan( 0,  2); //                                     [ 0,  2]
      sp.fill_plan(18, 20); //                                     [ 3,  5]
      sp.fill_plan(21, 26); // TID,  6 disks  (3 mono + 3 stereo)  [ 6, 11]
      sp.fill_plan(27, 44); // TEC, 18 disks  (9 mono + 9 stereo)  [12, 29]
      sp.set_iterator_limits(2, 0);
    }
  }

/*
  ////// Example backward search setup for initialStep iteration (currently 'replaced' by seed duplicate merging)
  void SetupBackwardSearch_Iter0(IterationConfig& ic)
  {
    ic.m_backward_params = ic.m_params;
    ic.m_backward_search = true;
    ic.m_backward_drop_seed_hits = true;
    ic.m_backward_fit_min_hits   = 8;
    auto &spv = ic.m_steering_params;
    spv[TrackerInfo::Reg_Endcap_Neg]    .set_iterator_limits(2, 3, 5);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(2, 3, 7);
    spv[TrackerInfo::Reg_Barrel]        .set_iterator_limits(2, 3, 4);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(2, 3, 7);
    spv[TrackerInfo::Reg_Endcap_Pos]    .set_iterator_limits(2, 3, 5);
  }
*/

  void SetupBackwardSearch_PixelCommon(IterationConfig& ic)
  {
    ic.m_backward_params = ic.m_params;
    ic.m_backward_search = true;
    ic.m_backward_drop_seed_hits = false;
    ic.m_backward_fit_min_hits   = 99;
    auto &spv = ic.m_steering_params;
    spv[TrackerInfo::Reg_Endcap_Neg]    .set_iterator_limits(2, 0, 3);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(2, 0, 4);
    spv[TrackerInfo::Reg_Barrel]        .set_iterator_limits(2, 0, 2);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(2, 0, 4);
    spv[TrackerInfo::Reg_Endcap_Pos]    .set_iterator_limits(2, 0, 3);
  }

  void SetupBackwardSearch_Iter7(IterationConfig& ic)
  {
    ic.m_backward_params = ic.m_params;
    ic.m_backward_search = true;
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
    ic.m_backward_params = ic.m_params;
    ic.m_backward_search = true;
    ic.m_backward_params.maxHolesPerCand = 2;
    ic.m_backward_params.maxConsecHoles  = 2;
    // Remove pixel/tib/tid layers from FwdSearch, add them to BkwSearch/
    auto &spv = ic.m_steering_params;
    spv[TrackerInfo::Reg_Endcap_Neg]    .set_iterator_limits(12, 12, 24);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(22, 19, 39);
    spv[TrackerInfo::Reg_Barrel]        .set_iterator_limits(12, 10, 14);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(22, 19, 39);
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
      ip.chi2Cut_min      = 15.0;
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    }
    else if (it == 1) // for triplet steps, nlayers_per_seed=3
    {   
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed  = 5;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 15.0;
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    }
    else if (it == 2)
    {   
      ip.nlayers_per_seed = 4;
      ip.maxCandsPerSeed  = 5;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 15.0;
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    } 
    else if (it == 3) // for triplet steps, nlayers_per_seed=3
    {   
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed  = 5;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 15.0;
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    } 
    else if (it == 4)
    {   
      ip.nlayers_per_seed = 4;
      ip.maxCandsPerSeed  = 5;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 15.0;
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    } 
    else if (it == 5) // for triplet steps, nlayers_per_seed=3
    {   
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed  = 5;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 15.0;
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    } 
    else if (it == 6) // for triplet steps, nlayers_per_seed=3; for mixeTripletSetp, also maxCandsPerSeed=2
    {
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed  = 2;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 15.0;
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    } 
    else if (it == 7) // for PixelLess step, maxCandsPerSeed=2 and maxHolesPerCand=maxConsecHoles=0
    {   
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed  = 2;
      ip.maxHolesPerCand  = 0;
      ip.maxConsecHoles   = 1;
      ip.chi2Cut_min      = 15.0;
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    }   
    else if (it == 8) // for TobTec step, maxCandsPerSeed=2 and maxHolesPerCand=maxConsecHoles=0
    {   
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed  = 2;
      ip.maxHolesPerCand  = 0;
      ip.maxConsecHoles   = 1;
      ip.chi2Cut_min      = 15.0;
      ip.chi2CutOverlap   = 3.5;
      ip.pTCutOverlap     = 1;
    }
    else if (it == 9) // addign also pixel pair step - algo -> 6
    {
      ip.nlayers_per_seed = 2;
      ip.maxCandsPerSeed  = 3;
      ip.maxHolesPerCand  = 4;
      ip.maxConsecHoles   = 2;
      ip.chi2Cut_min      = 15.0;
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


  //=================
  // partitionSeeds0
  //=================

  [[maybe_unused]] void partitionSeeds0(const TrackerInfo &trk_info,
                                        const TrackVec &in_seeds,
                                        const EventOfHits &eoh,
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
      part.m_sort_score[i] = 7.0f * (reg - 2) + eta;
    }
  }


  //=================
  // partitionSeeds1
  //=================

  [[maybe_unused]] void partitionSeeds1(const TrackerInfo &trk_info,
                                        const TrackVec &in_seeds,
                                        const EventOfHits &eoh,
                                        IterationSeedPartition &part) {
    // Seeds are placed into eta regions and sorted on region + eta.

    const LayerInfo &tib1 = trk_info.m_layers[4];
    //const LayerInfo &tib6 = trk_info.m_layers[9];
    const LayerInfo &tob1 = trk_info.m_layers[10];
    //const LayerInfo &tob8 = trk_info.m_layers[17];

    const LayerInfo &tidp1 = trk_info.m_layers[21];
    const LayerInfo &tidn1 = trk_info.m_layers[48];

    const LayerInfo &tecp1 = trk_info.m_layers[27];
    const LayerInfo &tecn1 = trk_info.m_layers[54];

    // Merge first two layers to account for mono/stereo coverage.
    // TrackerInfo could hold joint limits for sub-detectors.
    const auto &L = trk_info.m_layers;
    const float tidp_rin  = std::min(L[21].m_rin,  L[22].m_rin);
    const float tidp_rout = std::max(L[21].m_rout, L[22].m_rout);
    const float tecp_rin  = std::min(L[27].m_rin,  L[28].m_rin);
    const float tecp_rout = std::max(L[27].m_rout, L[28].m_rout);
    const float tidn_rin  = std::min(L[48].m_rin,  L[49].m_rin);
    const float tidn_rout = std::max(L[48].m_rout, L[49].m_rout);
    const float tecn_rin  = std::min(L[54].m_rin,  L[55].m_rin);
    const float tecn_rout = std::max(L[54].m_rout, L[55].m_rout);

    const float tid_z_extra = 0.0f; //  5.0f;
    const float tec_z_extra = 0.0f; // 10.0f;

    const int size = in_seeds.size();

    auto barrel_pos_check = [](const Track &S, float maxR, float rin, float zmax)->bool
    {
      bool inside = maxR > rin && S.zAtR(rin) < zmax;
      return inside;
    };

    auto barrel_neg_check = [](const Track &S, float maxR, float rin, float zmin)->bool
    {
      bool inside = maxR > rin && S.zAtR(rin) > zmin;
      return inside;
    };

    auto endcap_pos_check = [](const Track &S, float maxR, float rout, float rin, float zmin)->bool
    {
      bool inside = maxR > rout ?  S.zAtR(rout) > zmin :
                    (maxR > rin  && S.zAtR(maxR) > zmin);
      return inside;
    };

    auto endcap_neg_check = [](const Track &S, float maxR, float rout, float rin, float zmax)->bool
    {
      bool inside = maxR > rout ?  S.zAtR(rout) < zmax :
                    (maxR > rin  && S.zAtR(maxR) < zmax);
      return inside;
    };

    for (int i = 0; i < size; ++i) {
      const Track &S = in_seeds[i];

      HitOnTrack hot = S.getLastHitOnTrack();
      float eta = eoh[hot.layer].GetHit(hot.index).eta();
      // float  eta = S.momEta();

      // Region to be defined by propagation / intersection tests
      TrackerInfo::EtaRegion reg;

      const bool  z_dir_pos = S.pz() > 0;
      const float maxR = S.maxReachRadius();

      if (z_dir_pos)
      {
        bool in_tib = barrel_pos_check(S, maxR, tib1.m_rin, tib1.m_zmax);
        bool in_tob = barrel_pos_check(S, maxR, tob1.m_rin, tob1.m_zmax);

        if (!in_tib && !in_tob) {
          reg = TrackerInfo::Reg_Endcap_Pos;
        } else {
          bool in_tid = endcap_pos_check(S, maxR, tidp_rout, tidp_rin, tidp1.m_zmin - tid_z_extra);
          bool in_tec = endcap_pos_check(S, maxR, tecp_rout, tecp_rin, tecp1.m_zmin - tec_z_extra);

          if (!in_tid && !in_tec) {
            reg = TrackerInfo::Reg_Barrel;
          } else {
            reg = TrackerInfo::Reg_Transition_Pos;
          }
        }
      }
      else
      {
        bool in_tib = barrel_neg_check(S, maxR, tib1.m_rin, tib1.m_zmin);
        bool in_tob = barrel_neg_check(S, maxR, tob1.m_rin, tob1.m_zmin);

        if (!in_tib && !in_tob) {
          reg = TrackerInfo::Reg_Endcap_Neg;
        } else {
          bool in_tid = endcap_neg_check(S, maxR, tidn_rout, tidn_rin, tidn1.m_zmax + tid_z_extra);
          bool in_tec = endcap_neg_check(S, maxR, tecn_rout, tecn_rin, tecn1.m_zmax + tec_z_extra);

          if (!in_tid && !in_tec) {
            reg = TrackerInfo::Reg_Barrel;
          } else {
            reg = TrackerInfo::Reg_Transition_Neg;
          }
        }
      }

      part.m_region[i] = reg;
      part.m_sort_score[i] = 7.0f * (reg - 2) + eta;
    }
  }

  //======================
  // partitionSeeds1debug
  //======================

  [[maybe_unused]] void partitionSeeds1debug(const TrackerInfo &trk_info,
                                             const TrackVec &in_seeds,
                                             const EventOfHits &eoh,
                                             IterationSeedPartition &part)
  {
    // Seeds are placed into eta regions and sorted on region + eta.

    const LayerInfo &tib1 = trk_info.m_layers[4];
    //const LayerInfo &tib6 = trk_info.m_layers[9];
    const LayerInfo &tob1 = trk_info.m_layers[10];
    //const LayerInfo &tob8 = trk_info.m_layers[17];

    const LayerInfo &tidp1 = trk_info.m_layers[21];
    const LayerInfo &tidn1 = trk_info.m_layers[48];

    const LayerInfo &tecp1 = trk_info.m_layers[27];
    const LayerInfo &tecn1 = trk_info.m_layers[54];

    // Merge first two layers to account for mono/stereo coverage.
    // TrackerInfo could hold joint limits for sub-detectors.
    const auto &L = trk_info.m_layers;
    const float tidp_rin  = std::min(L[21].m_rin,  L[22].m_rin);
    const float tidp_rout = std::max(L[21].m_rout, L[22].m_rout);
    const float tecp_rin  = std::min(L[27].m_rin,  L[28].m_rin);
    const float tecp_rout = std::max(L[27].m_rout, L[28].m_rout);
    const float tidn_rin  = std::min(L[48].m_rin,  L[49].m_rin);
    const float tidn_rout = std::max(L[48].m_rout, L[49].m_rout);
    const float tecn_rin  = std::min(L[54].m_rin,  L[55].m_rin);
    const float tecn_rout = std::max(L[54].m_rout, L[55].m_rout);

    const float tid_z_extra = 0.0f; //  5.0f;
    const float tec_z_extra = 0.0f; // 10.0f;

    const int size = in_seeds.size();

    auto barrel_pos_check = [](const Track &S, float maxR, float rin, float zmax, const char *det)->bool
    {
      bool inside = maxR > rin && S.zAtR(rin) < zmax;

      printf("  in_%s=%d  maxR=%7.3f, rin=%7.3f -- ", det, inside, maxR, rin);
      if (maxR > rin) {
        printf("maxR > rin:   S.zAtR(rin) < zmax  -- %.3f <? %.3f\n", S.zAtR(rin), zmax);
      } else {
        printf("maxR < rin: no pie.\n");
      }

      return inside;
    };

    auto barrel_neg_check = [](const Track &S, float maxR, float rin, float zmin, const char *det)->bool
    {
      bool inside = maxR > rin && S.zAtR(rin) > zmin;

      printf("  in_%s=%d  maxR=%7.3f, rin=%7.3f -- ", det, inside, maxR, rin);
      if (maxR > rin) {
        printf("maxR > rin:   S.zAtR(rin) > zmin  -- %.3f >? %.3f\n", S.zAtR(rin), zmin);
      } else {
        printf("maxR < rin: no pie.\n");
      }

      return inside;
    };

    auto endcap_pos_check = [](const Track &S, float maxR, float rout, float rin, float zmin, const char *det)->bool
    {
      bool inside = maxR > rout ?  S.zAtR(rout) > zmin :
                    (maxR > rin  && S.zAtR(maxR) > zmin);

      printf("  in_%s=%d  maxR=%7.3f, rout=%7.3f, rin=%7.3f -- ", det, inside, maxR, rout, rin);
      if (maxR > rout) {
        printf("maxR > rout:  S.zAtR(rout) > zmin  -- %.3f >? %.3f\n", S.zAtR(rout), zmin);
      } else if (maxR > rin) {
        printf("maxR > rin:   S.zAtR(maxR) > zmin) -- %.3f >? %.3f\n", S.zAtR(maxR), zmin);
      } else {
        printf("maxR < rin: no pie.\n");
      }

      return inside;
    };

    auto endcap_neg_check = [](const Track &S, float maxR, float rout, float rin, float zmax, const char *det)->bool
    {
      bool inside = maxR > rout ?  S.zAtR(rout) < zmax :
                    (maxR > rin  && S.zAtR(maxR) < zmax);

      printf("  in_%s=%d  maxR=%7.3f, rout=%7.3f, rin=%7.3f -- ", det, inside, maxR, rout, rin);
      if (maxR > rout) {
        printf("maxR > rout:  S.zAtR(rout) < zmax  -- %.3f <? %.3f\n", S.zAtR(rout), zmax);
      } else if (maxR > rin) {
        printf("maxR > rin:   S.zAtR(maxR) < zmax  -- %.3f <? %.3f\n", S.zAtR(maxR), zmax);
      } else {
        printf("maxR < rin: no pie.\n");
      }

      return inside;
    };

    for (int i = 0; i < size; ++i) {
      const Track &S = in_seeds[i];

      HitOnTrack hot = S.getLastHitOnTrack();
      float eta = eoh[hot.layer].GetHit(hot.index).eta();
      // float  eta = S.momEta();

      // Region to be defined by propagation / intersection tests
      TrackerInfo::EtaRegion reg;

      const bool  z_dir_pos = S.pz() > 0;
      const float maxR = S.maxReachRadius();

      printf("partitionSeeds1debug seed index %d, z_dir_pos=%d (pz=%.3f), maxR=%.3f\n",
              i, z_dir_pos, S.pz(), maxR);

      if (z_dir_pos)
      {
        bool in_tib = barrel_pos_check(S, maxR, tib1.m_rin, tib1.m_zmax, "TIBp");
        bool in_tob = barrel_pos_check(S, maxR, tob1.m_rin, tob1.m_zmax, "TOBp");

        if (!in_tib && !in_tob) {
          reg = TrackerInfo::Reg_Endcap_Pos;
          printf("  --> region = %d, endcap pos\n", reg);
        } else {
          bool in_tid = endcap_pos_check(S, maxR, tidp_rout, tidp_rin, tidp1.m_zmin - tid_z_extra, "TIDp");
          bool in_tec = endcap_pos_check(S, maxR, tecp_rout, tecp_rin, tecp1.m_zmin - tec_z_extra, "TECp");

          if (!in_tid && !in_tec) {
            reg = TrackerInfo::Reg_Barrel;
            printf("  --> region = %d, barrel\n", reg);
          } else {
            reg = TrackerInfo::Reg_Transition_Pos;
            printf("  --> region = %d, transition pos\n", reg);
          }
        }
      }
      else
      {
        bool in_tib = barrel_neg_check(S, maxR, tib1.m_rin, tib1.m_zmin, "TIBn");
        bool in_tob = barrel_neg_check(S, maxR, tob1.m_rin, tob1.m_zmin, "TOBn");

        if (!in_tib && !in_tob) {
          reg = TrackerInfo::Reg_Endcap_Neg;
          printf("  --> region = %d, endcap neg\n", reg);
        } else {
          bool in_tid = endcap_neg_check(S, maxR, tidn_rout, tidn_rin, tidn1.m_zmax + tid_z_extra, "TIDn");
          bool in_tec = endcap_neg_check(S, maxR, tecn_rout, tecn_rin, tecn1.m_zmax + tec_z_extra, "TECn");

          if (!in_tid && !in_tec) {
            reg = TrackerInfo::Reg_Barrel;
            printf("  --> region = %d, barrel\n", reg);
          } else {
            reg = TrackerInfo::Reg_Transition_Neg;
            printf("  --> region = %d, transition neg\n", reg);
          }
        }
      }

      part.m_region[i] = reg;
      part.m_sort_score[i] = 7.0f * (reg - 2) + eta;
    }
  }


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
    ii[0].m_partition_seeds = partitionSeeds1;

    SetupCoreSteeringParams_Iter0(ii[0]);

    // At this point copy out layer/steering stuff for reuse in later iterations.
    IterationConfig def_itconf;
    def_itconf.CloneLayerSteerCore(ii[0]);

    SetupIterationParams(ii[0].m_params, 0);
    ii[0].set_dupclean_flag();
    ii[0].set_dupl_params(0.24, 0.002,0.004,0.008);
    fill_hit_selection_windows_params(ii[0]);
    ii[0].m_backward_params = ii[0].m_params;

    ii[1].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[1].m_params, 1);
    ii[1].set_iteration_index_and_track_algorithm(1, (int) TrackBase::TrackAlgorithm::highPtTripletStep);
    ii[1].set_seed_cleaning_params(2.0, 0.018, 0.018, 0.018, 0.018, 0.018, 0.05, 0.018, 0.05); 
    ii[1].set_dupclean_flag();
    ii[1].set_dupl_params(0.24, 0.03,0.05,0.08);
    fill_hit_selection_windows_params(ii[1]);
    SetupBackwardSearch_PixelCommon(ii[1]);

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
    SetupBackwardSearch_PixelCommon(ii[3]);
    
    ii[4].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[4].m_params, 4);
    ii[4].set_iteration_index_and_track_algorithm(4, (int) TrackBase::TrackAlgorithm::detachedQuadStep);
    ii[4].set_seed_cleaning_params(2.0, 0.018, 0.018, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05);
    ii[4].set_dupclean_flag();
    ii[4].set_dupl_params(0.24, 0.018,0.05,0.05);
    fill_hit_selection_windows_params(ii[4]);
    ii[4].m_backward_params = ii[4].m_params;
    
    ii[5].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[5].m_params, 5);
    ii[5].set_iteration_index_and_track_algorithm(5, (int) TrackBase::TrackAlgorithm::detachedTripletStep);
    ii[5].set_seed_cleaning_params(2.0, 0.018, 0.018, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05);
    ii[5].set_dupclean_flag();
    ii[5].set_dupl_params(0.24, 0.01,0.01,0.1);
    ii[5].m_requires_quality_filter = true;
    fill_hit_selection_windows_params(ii[5]);
    SetupBackwardSearch_PixelCommon(ii[5]);

    ii[6].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[6].m_params, 6);
    ii[6].set_iteration_index_and_track_algorithm(6, (int) TrackBase::TrackAlgorithm::mixedTripletStep);
    ii[6].set_seed_cleaning_params(2.0, 0.05, 0.05, 0.135, 0.135, 0.05, 0.05, 0.135, 0.135);
    ii[6].set_dupclean_flag();
    ii[6].set_dupl_params(0.2, 0.05,0.05,0.05); 
    fill_hit_selection_windows_params(ii[6]);
    SetupBackwardSearch_PixelCommon(ii[6]);

    ii[7].CloneLayerSteerCore(def_itconf);
    SetupIterationParams(ii[7].m_params, 7);
    ii[7].set_iteration_index_and_track_algorithm(7, (int) TrackBase::TrackAlgorithm::pixelLessStep);
    ii[7].set_seed_cleaning_params(2.0, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135);
    ii[7].set_qf_flags();
    ii[7].set_qf_params(3,0.14);
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
    ii[9].m_requires_quality_filter = true;
    fill_hit_selection_windows_params(ii[9]);
    SetupBackwardSearch_PixelCommon(ii[9]);

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
