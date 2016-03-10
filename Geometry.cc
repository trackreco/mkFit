#include "Geometry.h"

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

Geometry Geometry::clone() const
{
  Geometry NewGeo(*this);
  for (auto&& t : NewGeo.solids_) {
    t = t->Clone();
  }
  return NewGeo;
}
