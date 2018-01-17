#include "TrackerInfo.h"

#include <cassert>

void LayerInfo::set_limits(float r1, float r2, float z1, float z2)
{
  m_rin = r1; m_rout = r2; m_zmin = z1; m_zmax = z2;
}

void LayerInfo::set_next_layers(int nb, int nep, int nen)
{
  m_next_barrel = nb; m_next_ecap_pos = nep; m_next_ecap_neg = nen;
}

void LayerInfo::set_selection_limits(float p1, float p2, float q1, float q2)
{
  m_select_min_dphi = p1; m_select_max_dphi = p2;
  m_select_min_dq   = q1; m_select_max_dq   = q2;
}

void LayerInfo::set_r_hole_range(float rh1, float rh2)
{
  m_has_r_range_hole = true;
  m_hole_r_min = rh1;
  m_hole_r_max = rh2;
}


//==============================================================================
// TrackerInfo
//==============================================================================

LayerInfo TrackerInfo::s_undefined_layer(-1, LayerInfo::Undef);

void TrackerInfo::set_eta_regions(float tr_beg, float tr_end, float ec_end, bool has_sibl_lyrs)
{
  m_eta_trans_beg = tr_beg;
  m_eta_trans_end = tr_end;
  m_eta_ecap_end  = ec_end;
  m_has_sibling_layers = has_sibl_lyrs;
}

void TrackerInfo::reserve_layers(int n_brl, int n_ec_pos, int n_ec_neg)
{
  m_layers.reserve(n_brl + n_ec_pos + n_ec_neg);
  m_barrel.reserve(n_brl);
  m_ecap_pos.reserve(n_ec_pos);
  m_ecap_neg.reserve(n_ec_neg);
}

void TrackerInfo::create_layers(int n_brl, int n_ec_pos, int n_ec_neg)
{
  reserve_layers(n_brl, n_ec_pos, n_ec_neg);
  for (int i = 0; i < n_brl;    ++i) new_barrel_layer();
  for (int i = 0; i < n_ec_pos; ++i) new_ecap_pos_layer();
  for (int i = 0; i < n_ec_neg; ++i) new_ecap_neg_layer();
}

int TrackerInfo::new_layer(LayerInfo::LayerType_e type)
{
  int l = (int) m_layers.size();
  m_layers.emplace_back(LayerInfo(l, type));
  return l;
}

LayerInfo& TrackerInfo::new_barrel_layer()
{
  m_barrel.push_back( new_layer(LayerInfo::Barrel) );
  return m_layers.back();
}

LayerInfo& TrackerInfo::new_ecap_pos_layer()
{
  m_ecap_pos.push_back( new_layer(LayerInfo::EndCapPos) );
  return m_layers.back();
}

LayerInfo& TrackerInfo::new_ecap_neg_layer()
{
  m_ecap_neg.push_back( new_layer(LayerInfo::EndCapNeg) );
  return m_layers.back();
}

//------------------------------------------------------------------------------

bool TrackerInfo::are_layers_siblings(int l1, int l2) const
{
  assert(l1 < m_layers.size() && l2 < m_layers.size());

  const LayerInfo &i1 = m_layers[l1];
  const LayerInfo &i2 = m_layers[l2];

  if ( ! m_has_sibling_layers || i1.m_layer_type == i2.m_layer_type) return false;

  if (i1.m_layer_type == LayerInfo::Barrel)
    return l2 == i1.m_sibl_ecap_pos || l2 == i1.m_sibl_ecap_neg;
  else
    return l2 == i1.m_sibl_barrel;
}


//==============================================================================
// Plugin Loader
//==============================================================================

#include <dlfcn.h>
#include <sys/stat.h>
#include <cstdlib>

namespace
{
  const char *search_path[] = { "", "../Geoms/", "Geoms/", "../", 0 };
}

void TrackerInfo::ExecTrackerInfoCreatorPlugin(const std::string& base, TrackerInfo &ti, bool verbose)
{
#ifdef __MIC__
  std::string soname = base + "-mic.so";
#else
  std::string soname = base + ".so";
#endif

  struct stat st;

  int si = 0;
  while (search_path[si])
  {
    std::string path(search_path[si]); path += soname;
    if (stat(path.c_str(), &st) == 0)
    {
      printf("TrackerInfo::ExecTrackerInfoCreatorPlugin processing '%s'\n", path.c_str());

      void *h = dlopen(path.c_str(), RTLD_LAZY);
      if (!h) { perror("dlopen failed"); exit(2); }

      long long* p2f = (long long*) dlsym(h, "TrackerInfoCrator_ptr");
      if (!p2f) { perror("dlsym failed"); exit(2); }

      TrackerInfoCreator_foo foo = (TrackerInfoCreator_foo)(*p2f);
      foo(ti, verbose);

      // XXXXMT4Dan
      // With dlclose I saw (rarely) valgrind errors saying
      // Invalid read of size 4
      //     at 0x4119F7: is_barrel (TrackerInfo.h:54)
      // Commenting it out for now.
      // dlclose(h);
      return;
    }

    ++si;
  }

  fprintf(stderr, "TrackerInfo plugin '%s' not found in search path.\n", soname.c_str());
  exit(2);
}
