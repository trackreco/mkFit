#ifndef TRACKERINFO_H
#define TRACKERINFO_H

#include "Matrix.h"

#include <string>
#include <vector>
#include <stdexcept>

#include <cmath>
#include <cstdio>

//==============================================================================

enum WithinSensitiveRegion_e
{
  WSR_Undef = -1, WSR_Inside = 0, WSR_Edge, WSR_Outside
};

struct WSR_Result
{
  // Could also store XHitSize count equivalent here : 16;
  WithinSensitiveRegion_e m_wsr    : 8;
  bool                    m_in_gap : 8;

  WSR_Result() : m_wsr(WSR_Undef), m_in_gap(false) {}

  WSR_Result(WithinSensitiveRegion_e wsr, bool in_gap) : m_wsr(wsr), m_in_gap(in_gap) {}
};

//==============================================================================

class LayerInfo
{
private:
  bool  is_in_r_hole_no_check(float r) const { return r > m_hole_r_min && r < m_hole_r_max; }

public:
  enum  LayerType_e { Undef = -1, Barrel = 0, EndCapPos = 1, EndCapNeg = 2 };

  int           m_layer_id   = -1;
  LayerType_e   m_layer_type = Undef;

  float         m_rin, m_rout, m_zmin, m_zmax;
  float         m_propagate_to;

  int           m_next_barrel = -1, m_next_ecap_pos = -1, m_next_ecap_neg = -1;
  int           m_sibl_barrel = -1, m_sibl_ecap_pos = -1, m_sibl_ecap_neg = -1;

  bool          m_is_outer         = false;
  bool          m_has_r_range_hole = false;
  float         m_hole_r_min, m_hole_r_max; // This could be turned into std::function when needed.

  // Selection limits
  float         m_q_bin; // > 0 - bin width, < 0 - number of bins
  float         m_select_min_dphi, m_select_max_dphi;
  float         m_select_min_dq,   m_select_max_dq;

  // Additional stuff needed?
  // * pixel / strip, mono / stereo
  // * resolutions, min/max search windows
  // * holes in coverage
  // * functions / lambdas for deciding / calculating stuff
  // * ...
  // * pointers to hit containers

  LayerInfo(int lid, LayerType_e type) :
    m_layer_id(lid),
    m_layer_type(type)
  {}

  void  set_limits(float r1, float r2, float z1, float z2);
  void  set_next_layers(int nb, int nep, int nen);
  void  set_selection_limits(float p1, float p2, float q1, float q2);
  void  set_r_hole_range(float rh1, float rh2);

  float r_mean()    const { return (m_rin  + m_rout) / 2; }
  float z_mean()    const { return (m_zmin + m_zmax) / 2; }

  bool  is_barrel() const { return m_layer_type == Barrel; }

  bool  is_within_z_limits(float z) const { return z > m_zmin && z < m_zmax; }
  bool  is_within_r_limits(float r) const { return r > m_rin  && r < m_rout; }
  bool  is_within_q_limits(float q) const { return is_barrel() ? is_within_z_limits(q) : is_within_r_limits(q); }

  bool  is_in_r_hole      (float r) const { return m_has_r_range_hole ? is_in_r_hole_no_check(r) : false; }

  WSR_Result is_within_z_sensitive_region(float z, float dz) const
  {
    if (z > m_zmax + dz || z < m_zmin - dz)  return WSR_Result(WSR_Outside, false);
    if (z < m_zmax - dz && z > m_zmin + dz)  return WSR_Result(WSR_Inside,  false);
    return WSR_Result(WSR_Edge, false);
  }

  WSR_Result is_within_r_sensitive_region(float r, float dr) const
  {
    if (r > m_rout + dr || r < m_rin - dr)  return WSR_Result(WSR_Outside, false);
    if (r < m_rout - dr && r > m_rin + dr)
    {
      if (m_has_r_range_hole)
      {
        if (r < m_hole_r_max - dr && r > m_hole_r_min + dr)  return WSR_Result(WSR_Outside, true);
        if (r < m_hole_r_max + dr && r > m_hole_r_min - dr ) return WSR_Result(WSR_Edge,    true);
      }
      return WSR_Result(WSR_Inside, false);
    }
    return WSR_Result(WSR_Edge, false);
  }

  void  print_layer()
  {
     printf("Layer %2d  r(%7.4f, %7.4f) z(% 9.4f, % 9.4f) next(%2d, %2d, %2d) is_brl=%d is_outer=%d\n",
            m_layer_id, m_rin, m_rout, m_zmin, m_zmax,
            m_next_barrel, m_next_ecap_pos, m_next_ecap_neg,
            is_barrel(), m_is_outer);
  }
};

//==============================================================================

class TrackerInfo
{
private:
  int new_layer(LayerInfo::LayerType_e type);

public:
  enum AbsEtaRegion_e { AbsReg_Outside = -1, AbsReg_Barrel = 0, AbsReg_Transition = 1, AbsReg_Endcap = 2 };

  enum EtaRegion { Reg_Begin = 0, Reg_Endcap_Neg = 0, Reg_Transition_Neg, Reg_Barrel,
                   Reg_Transition_Pos, Reg_Endcap_Pos, Reg_End, Reg_Count = Reg_End };

  std::vector<LayerInfo> m_layers;
  static LayerInfo       s_undefined_layer;

  std::vector<int>       m_barrel;
  std::vector<int>       m_ecap_pos;
  std::vector<int>       m_ecap_neg;

  float  m_eta_trans_beg, m_eta_trans_end, m_eta_ecap_end;
  bool   m_has_sibling_layers;


  void        set_eta_regions(float tr_beg, float tr_end, float ec_end,
                              bool has_sibl_lyrs);
  void        reserve_layers(int n_brl, int n_ec_pos, int n_ec_neg);
  void        create_layers (int n_brl, int n_ec_pos, int n_ec_neg);
  LayerInfo & new_barrel_layer();
  LayerInfo & new_ecap_pos_layer();
  LayerInfo & new_ecap_neg_layer();

  bool are_layers_siblings(int l1, int l2) const;

  bool is_barrel(float eta) const
  {
    return std::abs(eta) < m_eta_trans_beg;
  }

  bool is_transition(float eta, float safety = 0) const
  {
    return std::abs(eta) >= m_eta_trans_beg - safety && std::abs(eta) <= m_eta_trans_end + safety;
  }

  bool is_endcap(float eta) const
  {
    return std::abs(eta) > m_eta_trans_end;
  }

  EtaRegion find_eta_region(float eta) const
  {
    if      (eta < -m_eta_trans_end) return Reg_Endcap_Neg;
    else if (eta < -m_eta_trans_beg) return Reg_Transition_Neg;
    else if (eta <  m_eta_trans_beg) return Reg_Barrel;
    else if (eta <  m_eta_trans_end) return Reg_Transition_Pos;
    else                             return Reg_Endcap_Pos;
  }

  EtaRegion find_region_of_layer(int l) const
  {
    // Assumes layers increase monotonically for barrel / encap.
    // Never returns Transition region.

    if (l <= m_barrel.back())   return Reg_Barrel;
    if (l <= m_ecap_pos.back()) return Reg_Endcap_Pos;
    return Reg_Endcap_Neg;
  }

  const LayerInfo& outer_barrel_layer()
  {
    return m_layers[m_barrel.back()];
  }

  const LayerInfo& next_barrel_layer(int layer)
  {
    int nb = m_layers[layer].m_next_barrel;
    return (nb >= 0) ? m_layers[nb] : s_undefined_layer;
  }

  static void ExecTrackerInfoCreatorPlugin(const std::string& base, TrackerInfo &ti, bool verbose=false);
};

typedef void (*TrackerInfoCreator_foo)(TrackerInfo&, bool verbose);

#endif
