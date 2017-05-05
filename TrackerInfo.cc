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

//==============================================================================

bool TrackerInfo::are_layers_siblings(int l1, int l2) const
{
  assert(l1 < m_layers.size() && l2 < m_layers.size());

  const LayerInfo &i1 = m_layers[l1];

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
  std::string soname = base + ".so";

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

      dlclose(h);
      return;
    }

    ++si;
  }

  fprintf(stderr, "TrackerInfo plugin '%s' not found in search path.\n", soname.c_str());
  exit(2);
}
