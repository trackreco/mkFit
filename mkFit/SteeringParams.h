#ifndef SteeringParams_h
#define SteeringParams_h

#include "Track.h"
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
  // Up-link to IterationConfig, so that IterationLayerConfig only can be passed afterwards
  // NOTE: this solution assumes that IterationLayerConfig is not shared among different iterations.
  const IterationConfig &m_iter_config;

  // A pointer to LayerOfHits, so that there is no need to pass this along together with IterationLayerConfig.
  const LayerOfHits* m_layer_of_hits;

  // Stuff to be moved out of LayerInfo:
  // Selection limits are moved out of LayerInfo, since may be iteration-specific.
  // E.g., for low-pT iterations selection limits may need to be looser (?)
  float         m_select_min_dphi, m_select_max_dphi;
  float         m_select_min_dq,   m_select_max_dq;

  void set_selection_limits(float p1, float p2, float q1, float q2)
  {
    m_select_min_dphi = p1; m_select_max_dphi = p2;
    m_select_min_dq   = q1; m_select_max_dq   = q2;
  }

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
// IterationParams
//==============================================================================

class IterationParams
{
  // Iteration-specific parameters are moved out of Config, and re-initialized in Geoms/CMS-2017.cc:
  int nlayers_per_seed  = 3;
  int maxCandsPerSeed   = 5;
  int maxHolesPerCand   = 4;
  int maxConsecHoles    = 1;
  float chi2Cut        = 30;
  float chi2CutOverlap; // default: 5; cmssw: 3.5
  float pTCutOverlap; // default: 0; cmssw: 1

  // NOTE: iteration params could actually become layer-dependent, e.g., chi2Cut could be larger for first layers (?)
};


//==============================================================================
// IterationConfig
//==============================================================================

class Event;

struct MkSeedPacket
{
  int m_seedEtaSeparators_[5];
  int m_seedMinLastLayer_[5];
  int m_seedMaxLastLayer_[5];
  std::vector<HitVec> m_layerHits_;
  TrackVec m_inseeds_;
  TrackVec m_outtrks_;

  //----------------------------------------------------------------------------

  MkSeedPacket(Event *evt) :
    m_seedEtaSeparators(evt->seedEtaSeparators_), m_seedMinLastLayer_(evt->seedMinLastLayer_), m_seedMaxLastLayer_(evt->seedMaxLastLayer_), m_layerHits_(evt->layerHits_), m_inseeds_(evt->seedTracks_) {}

};

class IterationConfig
{
public:

  //Up-link to event, to be used by import_seeds
  const Event *ev;
  
  // Iteration index:
  const unsigned int    m_iter;

  // TrackerInfo reference (can be shared among builders / iterations) is declared here so that it is easier to pass along:
  const TrackerInfo     &m_tracker_info;
  
  // Reference to iteration parameters:
  const IterationParams &m_params;
  
  std::vector<IterationLayerConfig> m_iter_layer_configs;
  
  // Steering params and regions are kept as they are, just moved out of MkBuilder.
  // Steering params and regions are iteration-specific, thus initialized in Geoms/CMS-2017.cc
  std::vector<SteeringParams> m_steering_params[5];
  std::vector<int>            m_regions;

  // Virtual function for 'import_seeds' (previously in MkBuilder):
  virtual void import_seeds(ev, const TrackerInfo &ti, const unsigned int it=0);
    
  //----------------------------------------------------------------------------
  IterationConfig(const TrackerInfo &ti, const IterationParams &ip, const unsigned int it=0) :
    m_iter(it), m_tracker_info(ti), m_params(ip) {}

};


//==============================================================================

// IterationConfig instance is created in Geoms/CMS-2017.cc, and is passed to MkBuilder constructor (as const reference).
// Internally, MkBuilder is passing 'int region' to its own functions: should be fine as is.
// When calling MkFinder functions, MkBuilder is now passing a reference to LayerOfHits;
// this is now replaced by a reference to IterationLayerConfig, which includes ref's to IterationConfig and LayerOfHits.
// (Alternatively, one could pass both ref's to IterationLayerConfig and IterationConfig)
// Pointers to propagation functions can be passed as they are now.

} // end namespace mkfit
#endif
