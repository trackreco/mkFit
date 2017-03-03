#ifndef TRACKERINFO_H
#define TRACKERINFO_H

#include <vector>

#include <cstdio>

class LayerInfo
{
public:
  int   m_layer_id;

  float m_rin, m_rout, m_zmin, m_zmax;

  int   m_next_barrel, m_next_ecap_pos, m_next_ecap_neg;

  bool  m_is_barrel; // or byte enum barrel, ecappos, ecapneg
  bool  m_is_outer;

  // Additional stuff needed?
  // * pixel / strip, mono / stereo
  // * resolutions
  // * functions / lambdas for deciding / calculating stuff
  // * ...
  // * pointers to hit containers

  LayerInfo(int lid) : m_layer_id(lid) {}

  float r_mean() const { return (m_rin  + m_rout) / 2; }
  float z_mean() const { return (m_zmin + m_zmax) / 2; }

  void print_layer()
  {
     printf("Layer %2d  r(%7.4f, %7.4f) z(% 9.4f, % 9.4f) next(%2d, %2d, %2d) is_brl=%d is_outer=%d\n",
            m_layer_id, m_rin, m_rout, m_zmin, m_zmax,
            m_next_barrel, m_next_ecap_pos, m_next_ecap_neg,
            m_is_barrel, m_is_outer);
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
  std::vector<LayerInfo> m_layers;

  std::vector<int>       m_barrel;
  std::vector<int>       m_ecap_pos;
  std::vector<int>       m_ecap_neg;

  float       m_eta_trans_beg, m_eta_trans_end, m_eta_max;

  LayerInfo & new_barrel_layer()   { m_barrel  .push_back( new_layer() ); return m_layers.back(); }
  LayerInfo & new_ecap_pos_layer() { m_ecap_pos.push_back( new_layer() ); return m_layers.back(); }
  LayerInfo & new_ecap_neg_layer() { m_ecap_neg.push_back( new_layer() ); return m_layers.back(); }
};

#endif
