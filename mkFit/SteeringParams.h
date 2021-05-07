#ifndef SteeringParams_h
#define SteeringParams_h

namespace mkfit {

//==============================================================================
// LayerControl
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

} // end namespace mkfit

#endif
