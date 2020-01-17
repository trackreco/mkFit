#ifndef SteeringParams_h
#define SteeringParams_h

#include "Matrix.h"

namespace mkfit {

class LayerInfo;
class CandCloner;
class MkBase;
class MkFitter;
class MkFinder;

#define COMPUTE_CHI2_ARGS const MPlexLS &,  const MPlexLV &, const MPlexQI &, \
                          const MPlexHS &,  const MPlexHV &, \
                          MPlexQF &,  const int, const PropagationFlags

#define UPDATE_PARAM_ARGS const MPlexLS &,  const MPlexLV &, MPlexQI &, \
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
// SteeringParams
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


//==============================================================================
// IterationLayerConfig
//==============================================================================

class IterationConfig;

class IterationLayerConfig
{
public:
  // Up-link to IterationConfig ... so we can pass this guy only.
  // Is this really smart? Can IterationLayerConfigs be shared among
  // several iterations? If yes, this will not work.
  // See also text at the end of file.
  const IterationConfig &m_iter_config;

  // Here so it does not have to passed along (is this reasonable?)
  // It might need to be re-set for every event (I forgot), that's why I have ptr,
  // otherwise ref would be preferable.
  const LayerOfHits* m_layer_of_hits;

  // Stuff to be moved out of LayerInfo:

  // Selection limits
  float         m_select_min_dphi, m_select_max_dphi;
  float         m_select_min_dq,   m_select_max_dq;

  // Adding hit selection limits dynamic factors
  float         m_qf_treg       = 1.0f;
  float         m_phif_treg     = 1.0f;
  float         m_phif_lpt_brl  = 1.0f;
  float         m_phif_lpt_treg = 1.0f;
  float         m_phif_lpt_ec   = 1.0f;

  //----------------------------------------------------------------------------

  IterationLayerConfig(const IterationConfig &ic) :
    m_iter_config (ic)
  {}
};


//==============================================================================
// IterationLayerConfig
//==============================================================================

class IterationParams
{
  // Stuff moved out from Config, like:

  int maxCandsPerSeed  = 6; // cmssw tests: 6 (GC had 3) \_ set from geom plugin
  int maxHolesPerCand  = 2; // cmssw tests: 12           /

  float chi2Cut = 15.;

  // Some iteration params could actually become layer-dependent, e.g.,
  // chi2Cut could be larger for first couple of layers.
};


//==============================================================================
// IterationConfig
//==============================================================================

class IterationConfig
{
public:
  // ptr or ref to TrackerInfo (can be shared among builders / iterations).
  // declared here so it's easier to pass along
  const TrackerInfo     &m_tracker_info;

  // ptr or ref to iteration parameters (can be shared)
  // Or just keep them here.
  const IterationParams &m_params;

  std::vector<IterationLayerConfig> m_iter_layer_configs;

  // LayerOfHits (in HitStructures.h) should NOT hold forwarding functions
  // into layer_info.

  // LayerOfHits should be passed in parallel ... or, be pointed to by IterationLayerConfig,

  // Steering params can be kept as they are, just moved out of MkBuilder.
  std::vector<SteeringParams> m_steering_params[5];
  std::vector<int>            m_regions;

  // It seems we are a little bit stuck with 5 regions for now.
  // We probably could loosen it up but, IIRC, there are some
  // checks in the code for region type, in particular for barrel.

  //----------------------------------------------------------------------------

  IterationConfig(const TrackerInfo &ti, const IterationParams &ip) :
    m_tracker_info (ti),
    m_params       (ip)
  {}

  // Here we also need either:
  // a) a virtual function that performs partitioning of seeds into regions; or
  // b) a std::function that takes MkBuilder, Event, and IterationConfig (or
  //    some subset of those objects / their components) and does the seed partitioning.
  //
  // Next question is how we want to do seed cleaning for multiple iterations.
};


//==============================================================================

// Then one would create IterationConfigs in the "geometry" plugin (or in
// CMSSW configuration processing class) and pass it to MkBuilder constructor,
// and MkBuilder would store it in a const ref.
//
// Internally MkBuilder passes 'int region' to it's own functions; this should still work.
//
// When calling MkFinder functions, MkBuilder now passes LayerOfHits& ... this will have to be
// changed to IterationLayerConfig& (if we include ref to IterationConfig in IterationLayerConfig)
// or to both IterationConfig& and IterationLayerConfig&.
// Pointers to propagation functions can be passed as they are now.

} // end namespace mkfit
#endif
