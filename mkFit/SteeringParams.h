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
                                MPlexQF &,  const int

#define UPDATE_PARAM_ARGS const MPlexLS &,  const MPlexLV &, const MPlexQI &, \
                          const MPlexHS &,  const MPlexHV &, \
                                MPlexLS &,        MPlexLV &, const int

struct FindingFoos
{
  void (*m_compute_chi2_foo)      (COMPUTE_CHI2_ARGS);
  void (*m_update_param_foo)      (UPDATE_PARAM_ARGS);
  void (MkBase::*m_propagate_foo) (float, const int);

  FindingFoos() {}

  FindingFoos(void (*cch2_f)      (COMPUTE_CHI2_ARGS),
              void (*updp_f)      (UPDATE_PARAM_ARGS),
              void (MkBase::*p_f) (float, const int)) :
    m_compute_chi2_foo(cch2_f),
    m_update_param_foo(updp_f),
    m_propagate_foo(p_f)
  {}

};

struct LayerControl
{
  int  m_layer;

  // Idea only ... need some parallel structure for candidates to make sense (where i can store it).
  // Or have per layer containers where I place track indices to enable. Or something. Sigh.
  int  m_on_miss_jump_to;
  int  m_on_hit_jump_to;

  bool m_pickup_only; // do not propagate to this layer and process hits, pickup seeds only.

  LayerControl(int lay, bool pu_only=false) : m_layer(lay), m_pickup_only(pu_only) {}
};

struct SteeringParams
{
  void (*compute_chi2_foo) (COMPUTE_CHI2_ARGS);
  void (*update_param_foo) (UPDATE_PARAM_ARGS);

  int   LayerInfo::*next_layer_doo;

  void (MkBase::*propagate_foo) (float, const int);

  std::vector<LayerControl> layer_plan;

  //----------------------------------------------------------------------------

  SteeringParams() {}

  SteeringParams(void (*cch2_f)(COMPUTE_CHI2_ARGS),
                 void (*updp_f)(UPDATE_PARAM_ARGS),
                 int   LayerInfo::*nl_d,
                 void (MkBase::*p_f)     (float, const int)) :
    compute_chi2_foo(cch2_f),
    update_param_foo(updp_f),
    next_layer_doo(nl_d),
    propagate_foo(p_f)
  {}

  void reserve_plan(int n)
  {
    layer_plan.reserve(n);
  }

  void append_plan(int layer, bool pu_only=false)
  {
    layer_plan.emplace_back(LayerControl(layer, pu_only));
  }

  void fill_plan(int first, int last, bool pu_only=false)
  {
    for (int i = first; i <= last; ++i)  append_plan(i, pu_only);
  }
};

#endif
