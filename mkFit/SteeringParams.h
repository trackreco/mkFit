#ifndef SteeringParams_h
#define SteeringParams_h

#include "Matrix.h"

#include "nlohmann/json_fwd.hpp"

#include <functional>

namespace mkfit {

class LayerInfo;
class CandCloner;
class MkBase;
class MkFitter;
class MkFinder;
class EventOfHits;

class Track;
typedef std::vector<Track>    TrackVec;

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
  // int  m_on_miss_jump_to = -999;
  // int  m_on_hit_jump_to  = -999;

  bool m_pickup_only;  // do not propagate to this layer and process hits, pickup seeds only.
  bool m_backfit_only; // layer only used in backward fit.

  //----------------------------------------------------------------------------

  LayerControl(int lay, bool pu_only=false, bool bf_only=false) :
    m_layer(lay), m_pickup_only(pu_only), m_backfit_only(bf_only) {}
};


//==============================================================================
// Hit masks / IterationMaskIfc
//==============================================================================

struct IterationMaskIfcBase
{
  virtual ~IterationMaskIfcBase() {}

  virtual const std::vector<bool>* get_mask_for_layer(int layer) const { return nullptr; }
};

struct IterationMaskIfc : public IterationMaskIfcBase
{
  std::vector<std::vector<bool>> m_mask_vector;

  const std::vector<bool>* get_mask_for_layer(int layer) const { return & m_mask_vector[layer]; }
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
  // Selection limits.
  float         m_select_min_dphi;
  float         m_select_max_dphi;
  float         m_select_min_dq;
  float         m_select_max_dq;

  void set_selection_limits(float p1, float p2, float q1, float q2)
  {
    m_select_min_dphi = p1; m_select_max_dphi = p2;
    m_select_min_dq   = q1; m_select_max_dq   = q2;
  }

  using dynamic_windows_foo = void(const Track &track, const float track_pt, const float track_eta,
                                   float &min_dq, float &max_dphi);

  std::function<dynamic_windows_foo> m_dynamic_windows;

  // Adding hit selection limits dynamic factors
  float         m_qf_treg       = 1.0f;
  float         m_phif_treg     = 1.0f;
  float         m_phif_lpt_brl  = 1.0f;
  float         m_phif_lpt_treg = 1.0f;
  float         m_phif_lpt_ec   = 1.0f;

  //----------------------------------------------------------------------------

  float min_dphi() const { return m_select_min_dphi; }
  float max_dphi() const { return m_select_max_dphi; }
  float min_dq()   const { return m_select_min_dq;   }
  float max_dq()   const { return m_select_max_dq;   }


  //----------------------------------------------------------------------------

  IterationLayerConfig() {}
};


//==============================================================================
// IterationParams
//==============================================================================

class IterationParams
{
public:
  // Iteration-specific parameters are moved out of Config, and re-initialized in Geoms/CMS-2017.cc:
  int   nlayers_per_seed  = 3;
  int   maxCandsPerSeed   = 5;
  int   maxHolesPerCand   = 4;
  int   maxConsecHoles    = 1;
  float chi2Cut           = 30;
  float chi2CutOverlap    = 3.5;
  float pTCutOverlap      = 1.0;
  // NOTE: iteration params could actually become layer-dependent, e.g., chi2Cut could be larger for first layers (?)
  
  //seed cleaning params
  float c_ptthr_hpt = 2.0;
  //initial
  float c_drmax_bh = 0.010;
  float c_dzmax_bh = 0.005;
  float c_drmax_eh = 0.010;
  float c_dzmax_eh = 0.010;
  float c_drmax_bl = 0.010;
  float c_dzmax_bl = 0.005;
  float c_drmax_el = 0.015;
  float c_dzmax_el = 0.015;
};


//==============================================================================
// IterationSeedPartition
//==============================================================================

class IterationSeedPartition
{
public:
  std::vector<int>    m_region;
  std::vector<float>  m_sort_score;

  IterationSeedPartition(int size) : m_region(size), m_sort_score(size) {}
};


//==============================================================================
// IterationConfig
//==============================================================================

class IterationConfig
{
public:
  using partition_seeds_foo = void(const TrackerInfo &, const TrackVec &, const EventOfHits &,
                                   IterationSeedPartition &);

  int    m_iteration_index  = -1;
  int    m_track_algorithm  = -1;

  // Iteration parameters (could be a ptr)
  IterationParams                     m_params;

  int                                 m_n_regions;
  std::vector<int>                    m_region_order;
  std::vector<SteeringParams>         m_steering_params;
  std::vector<IterationLayerConfig>   m_layer_configs;

  std::function<partition_seeds_foo>  m_partition_seeds;

  //----------------------------------------------------------------------------

  IterationConfig() {}

  void Clone(const IterationConfig &o)
  {
      // Clone iteration. m_iteration_index and m_track_algorithm are not copied
      // and need to be set separately.

      m_params          = o.m_params;

      m_n_regions       = o.m_n_regions;
      m_region_order    = o.m_region_order;
      m_steering_params = o.m_steering_params;
      m_layer_configs   = o.m_layer_configs;

      m_partition_seeds = o.m_partition_seeds;
  }

  void set_iteration_index_and_track_algorithm(int idx, int trk_alg)
  {
    m_iteration_index = idx;
    m_track_algorithm = trk_alg;
  }
  
  void set_seed_cleaning_params(float pt_thr, 
 				float dzmax_bh, float drmax_bh, 
				float dzmax_bl, float drmax_bl,
				float dzmax_eh, float drmax_eh,
				float dzmax_el, float drmax_el
				)
  {
       m_params.c_ptthr_hpt = pt_thr;
       m_params.c_drmax_bh = drmax_bh;
       m_params.c_dzmax_bh = dzmax_bh;
       m_params.c_drmax_eh = drmax_eh;
       m_params.c_dzmax_eh = dzmax_eh;
       m_params.c_drmax_bl = drmax_bl;
       m_params.c_dzmax_bl = dzmax_bl;
       m_params.c_drmax_el = drmax_el;
       m_params.c_dzmax_el = dzmax_el;
   }

  void set_num_regions_layers(int nreg, int nlay)
  {
    m_n_regions = nreg;
    m_region_order.resize(nreg);
    m_steering_params.resize(nreg);
    m_layer_configs.resize(nlay);
  }

  IterationLayerConfig& layer(int i) { return m_layer_configs[i]; }

  SteeringParams&       steering_params(int region) { return m_steering_params[region]; }
};


//==============================================================================
// IterationsInfo
//==============================================================================

class IterationsInfo
{
public:
  std::vector<IterationConfig> m_iterations;

  IterationsInfo() {}

  void resize(int ni) { m_iterations.resize(ni); }

  int  size() const { return m_iterations.size(); }

  IterationConfig& operator[](int i) { return m_iterations[i]; }
  const IterationConfig& operator[](int i) const { return m_iterations[i]; }
};


//==============================================================================

// IterationConfig instances are created in Geoms/CMS-2017.cc, Create_CMS_2017(),
// filling the IterationsInfo object passed in by reference.


//==============================================================================
// JSON config interface
//==============================================================================

class ConfigJsonPatcher
{
  nlohmann::json *m_json    = nullptr;
  nlohmann::json *m_current = nullptr;

  // add stack and cd_up() ? also, name stack for exceptions and printouts
  std::vector<nlohmann::json*> m_json_stack;
  std::vector<std::string>     m_path_stack;

  bool            m_owner   = false;
  bool            m_verbose = false;

  void release_json();

  std::string get_abs_path() const;
  std::string exc_hdr(const char* func=0) const;

public:
  ConfigJsonPatcher(bool verbose=false);
  ~ConfigJsonPatcher();

  template<class T> void Load(const T& o);
  template<class T> void Save(T& o);

  void cd    (const std::string& path);
  void cd_up (const std::string& path = "");
  void cd_top(const std::string& path = "");

  template<typename T>
  void replace(const std::string& path, T val);

  template<typename T>
  void replace(int first, int last, const std::string& path, T val);

  nlohmann::json& get(const std::string& path);

  int replace(const nlohmann::json &j);

  std::string dump(int indent=2);
};

void ConfigJson_Patch_File(IterationsInfo &its_info, const std::string &fname);

void ConfigJson_Test_Direct(IterationConfig &it_cfg);
void ConfigJson_Test_Patcher(IterationConfig &it_cfg);


} // end namespace mkfit

#endif
