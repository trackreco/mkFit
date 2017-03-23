#include "TrackerInfo.h"

#include <cassert>

bool TrackerInfo::are_layers_siblings(int l1, int l2) const
{
  assert(l1 < m_layers.size() && l2 < m_layers.size());

  const LayerInfo &i1 = m_layers[l1];

  if (i1.m_layer_type == LayerInfo::Barrel)
    return l2 == i1.m_sibl_ecap_pos || l2 == i1.m_sibl_ecap_neg;
  else
    return l2 == i1.m_sibl_barrel;
}
