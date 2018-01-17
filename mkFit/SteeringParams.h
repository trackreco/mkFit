#ifndef SteeringParams_h
#define SteeringParams_h

#include "Matrix.h"

class LayerInfo;
class CandCloner;
class MkBase;
class MkFitter;
class MkFinder;

#define COMPUTE_CHI2_ARGS const MPlexLS &,  const MPlexLV &, const MPlexQI &, \
                          const MPlexHS &,  const MPlexHV &, \
                          MPlexQF &,  const int, const PropagationFlags

#define UPDATE_PARAM_ARGS const MPlexLS &,  const MPlexLV &, const MPlexQI &, \
                          const MPlexHS &,  const MPlexHV &, \
                                MPlexLS &,        MPlexLV &, const int, const PropagationFlags

class FindingFoos
{
public:
  void (*m_compute_chi2_foo)      (COMPUTE_CHI2_ARGS);
  void (*m_update_param_foo)      (UPDATE_PARAM_ARGS);
  void (MkBase::*m_propagate_foo) (float, const int, const PropagationFlags);

  FindingFoos() {}

  FindingFoos(void (*cch2_f)      (COMPUTE_CHI2_ARGS),
              void (*updp_f)      (UPDATE_PARAM_ARGS),
              void (MkBase::*p_f) (float, const int, const PropagationFlags)) :
    m_compute_chi2_foo(cch2_f),
    m_update_param_foo(updp_f),
    m_propagate_foo(p_f)
  {}

};

//==============================================================================
//==============================================================================

struct LayerControl
{
  int  m_layer;

  // Idea only ... need some parallel structure for candidates to make sense (where i can store it).
  // Or have per layer containers where I place track indices to enable. Or something. Sigh.
  int  m_on_miss_jump_to;
  int  m_on_hit_jump_to;

  bool m_pickup_only;  // do not propagate to this layer and process hits, pickup seeds only.
  bool m_backfit_only; // layer only used in backward fit.

  //----------------------------------------------------------------------------

  LayerControl(int lay, bool pu_only=false, bool bf_only=false) :
    m_layer(lay), m_pickup_only(pu_only), m_backfit_only(bf_only) {}
};

//==============================================================================

class SteeringParams
{
public:
  std::vector<LayerControl>           m_layer_plan;
  std::vector<LayerControl>::iterator m_begin_for_finding;
  std::vector<LayerControl>::iterator m_end_for_finding;

  //----------------------------------------------------------------------------

  SteeringParams() {}

  void reserve_plan(int n)
  {
    m_layer_plan.reserve(n);
  }

  void append_plan(int layer, bool pu_only=false, bool bf_only=false)
  {
    m_layer_plan.emplace_back(LayerControl(layer, pu_only, bf_only));
  }

  void fill_plan(int first, int last, bool pu_only=false, bool bf_only=false)
  {
    for (int i = first; i <= last; ++i)  append_plan(i, pu_only, bf_only);
  }

  void finalize_plan()
  {
    m_begin_for_finding = m_layer_plan.begin();
    while (m_begin_for_finding->m_backfit_only) ++m_begin_for_finding;
    m_end_for_finding = m_layer_plan.end();
  }

  std::vector<LayerControl>::iterator finding_begin() const { return m_begin_for_finding; }
  std::vector<LayerControl>::iterator finding_end()   const { return m_end_for_finding; }
};

#endif
