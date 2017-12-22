#include "Geometry.h"
//#define DEBUG
#include "Debug.h"

VUSolid::EnumInside Geometry::Inside (const UVector3 &aPoint) const {
  VUSolid::EnumInside in = VUSolid::eOutside;
  for (auto i = solids_.begin(); i != solids_.end(); ++i) {
    VUSolid::EnumInside t = (*i)->Inside(aPoint);
    if (t < in) in = t;
  }
  return in;
}

const VUSolid* Geometry::InsideWhat(const UVector3 &aPoint) const {
  for (auto i = solids_.begin(); i != solids_.end(); ++i) {
    if ((*i)->Inside(aPoint) != VUSolid::eOutside) {
      return *i;
    }
  }
  return nullptr;
}

int Geometry::LayerIndex(const UVector3 &aPoint) const {
  for (auto i = solids_.begin(); i != solids_.end(); ++i) {
    if ((*i)->Inside(aPoint) != VUSolid::eOutside) {
      return std::distance(solids_.begin(),i);
    }
  }
  return -1;
}

int Geometry::LayerOfSolid(const VUSolid *s) const
{
  for (auto i = solids_.begin(); i != solids_.end(); ++i) {
    if (*i == s) {
      return std::distance(solids_.begin(),i);
    }
  }
  return -1;
}

double Geometry::SafetyFromInside ( const UVector3 &aPoint, bool aAccurate) const {
  double small = 1e10;

  for (auto i = solids_.begin(); i != solids_.end(); ++i) {
    if ((*i)->Inside(aPoint) != VUSolid::eOutside) {
      double next = (*i)->SafetyFromInside(aPoint, aAccurate);
      if (next < small) small = next;
    }
  }
  return small;
}

double Geometry::SafetyFromOutside ( const UVector3 &aPoint, bool aAccurate) const {
  double small = 1e10;

  for (auto i = solids_.begin(); i != solids_.end(); ++i) {
    if ((*i)->Inside(aPoint) == VUSolid::eOutside) {
      double next = (*i)->SafetyFromOutside(aPoint, aAccurate);
      if (next < small) small = next;
    }
  }
  return small;
}

double Geometry::SafetyFromOutsideDr(const UVector3 &aPoint, double ooaCtgTheta,
                                     int skip_layer, int &layer, bool aAccurate) const
{
  //bool debug = false;

  double small = 1e10;

  dprintf("Geometry::SafetyFromOutsideDr r=%f, z=%f\n", aPoint.Perp(), aPoint.Z());
  int ii = 0;
  layer = -1;
  for (auto i = solids_.begin(); i != solids_.end(); ++i, ++ii)
  {
    if (ii != skip_layer && (*i)->Inside(aPoint) == VUSolid::eOutside)
    {
      double next = (*i)->SafetyFromOutsideDr(aPoint, ooaCtgTheta, aAccurate);
#ifdef DEBUG
      if (next < 15)  dprintf("    Radial distance to %2d = %4.5f\n", ii, next);
#endif
      if (next < small) { small = next; layer = ii; }
    }
  }
  dprintf("  Selected layer       %2d = %f\n", layer, small);
  return small;
}

Geometry Geometry::clone() const
{
  Geometry NewGeo(*this);
  for (auto&& t : NewGeo.solids_) {
    t = t->Clone();
  }
  return NewGeo;
}

//==============================================================================

#include "TrackerInfo.h"

void Geometry::BuildFromTrackerInfo(const TrackerInfo& tracker_info)
{
#ifndef WITH_USOLIDS
  for (auto& li : tracker_info.m_layers)
  {
    VUSolid* utub = new VUSolid(li.m_rin,  li.m_rout,
                                li.m_zmin, li.m_zmax,
                                li.is_barrel(), li.m_is_outer);
    AddLayer(utub, li.r_mean(), li.z_mean());
  }
#else
  // XXMT4D What do we do here?
  fprintf(stderr, "Geometry::BuildFromTrackerInfo only supported for SimpleGeometry.\n");
  exit(1);
#endif
}
