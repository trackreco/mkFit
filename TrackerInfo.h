#ifndef TRACKERINFO_H
#define TRACKERINFO_H

#include <vector>

#include <cmath>
#include <cstdio>

class LayerInfo
{
public:
  enum  LayerType_e { Barrel = 0, EndCapPos = 1, EndCapNeg = 2 };

  int           m_layer_id;
  LayerType_e   m_layer_type;

  float         m_rin, m_rout, m_zmin, m_zmax;

  int           m_next_barrel, m_next_ecap_pos, m_next_ecap_neg;
  int           m_sibl_barrel, m_sibl_ecap_pos, m_sibl_ecap_neg;

  bool          m_is_outer;

  // Additional stuff needed?
  // * pixel / strip, mono / stereo
  // * resolutions
  // * functions / lambdas for deciding / calculating stuff
  // * ...
  // * pointers to hit containers

  LayerInfo(int lid) : m_layer_id(lid) {}

  float r_mean()    const { return (m_rin  + m_rout) / 2; }
  float z_mean()    const { return (m_zmin + m_zmax) / 2; }

  bool  is_barrel() const { return m_layer_type == Barrel; }

  void  print_layer()
  {
     printf("Layer %2d  r(%7.4f, %7.4f) z(% 9.4f, % 9.4f) next(%2d, %2d, %2d) is_brl=%d is_outer=%d\n",
            m_layer_id, m_rin, m_rout, m_zmin, m_zmax,
            m_next_barrel, m_next_ecap_pos, m_next_ecap_neg,
            is_barrel(), m_is_outer);
  }
};

class TrackerInfo
{
private:
  int new_layer()
  {
    int l = (int) m_layers.size();
    m_layers.push_back(LayerInfo(l));
    return l;
  }

public:
  enum AbsEtaRegion_e { AbsReg_Outside = -1, AbsReg_Barrel = 0, AbsReg_Transition = 1, AbsReg_Endcap = 2 };

  enum EtaRegion { Reg_Begin = 0, Reg_Endcap_Neg = 0, Reg_Transition_Neg, Reg_Barrel,
                   Reg_Transition_Pos, Reg_Endcap_Pos, Reg_End, Reg_Count = Reg_End };

  std::vector<LayerInfo> m_layers;

  std::vector<int>       m_barrel;
  std::vector<int>       m_ecap_pos;
  std::vector<int>       m_ecap_neg;

  float  m_eta_trans_beg, m_eta_trans_end, m_eta_ecap_end;

  void set_eta_regions(float tr_beg, float tr_end, float ec_end)
  { m_eta_trans_beg = tr_beg; m_eta_trans_end = tr_end; m_eta_ecap_end = ec_end; }

  LayerInfo & new_barrel_layer()   { m_barrel  .push_back( new_layer() ); return m_layers.back(); }
  LayerInfo & new_ecap_pos_layer() { m_ecap_pos.push_back( new_layer() ); return m_layers.back(); }
  LayerInfo & new_ecap_neg_layer() { m_ecap_neg.push_back( new_layer() ); return m_layers.back(); }

  bool are_layers_siblings(int l1, int l2) const;

  bool is_barrel(float eta) const
  {
    return std::abs(eta) < m_eta_trans_beg;
  }

  bool is_transition(float eta) const
  {
    return std::abs(eta) >= m_eta_trans_beg && std::abs(eta) <= m_eta_trans_end;
  }

  bool is_endcap(float eta) const
  {
    return std::abs(eta) > m_eta_trans_end;
  }

  EtaRegion find_eta_region(float eta) const
  {
    if      (eta < -m_eta_trans_beg) return Reg_Endcap_Neg;
    else if (eta < -m_eta_trans_end) return Reg_Transition_Neg;
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
};

#endif
