#ifndef _geometry_
#define _geometry_

#include <vector>
#include "Matrix.h"

#ifdef WITH_USOLIDS
#include "USolids/include/VUSolid.hh"
#else
#include "SimpleGeom.h"
#endif

class Geometry {
public:
  void AddLayer(const VUSolid* s) { solids_.push_back(s); }
  void AddLayer(const VUSolid* s, const float radius, const float zPlane) { 
    solids_.push_back(s); radii_.push_back(radius); zPlanes_.push_back(zPlane);
  }
  unsigned int CountLayers() const { return solids_.size(); }
  const VUSolid* Layer(unsigned int l) const { return solids_.at(l); }
  float Radius(unsigned int l) const { return radii_.at(l); }
  float zPlane(unsigned int l) const { return zPlanes_.at(l); }
  VUSolid::EnumInside Inside (const UVector3 &aPoint) const;
  const VUSolid* InsideWhat(const UVector3 &aPoint) const;
  double  SafetyFromInside ( const UVector3 &aPoint, bool aAccurate=false) const;
  double  SafetyFromOutside( const UVector3 &aPoint, bool aAccurate=false) const;
private:
  std::vector<const VUSolid*> solids_;
  std::vector<float> radii_;
  std::vector<float> zPlanes_;
};
#endif
