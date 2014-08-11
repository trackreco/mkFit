#ifndef _geometry_
#define _geometry_

#include <vector>
#include "USolids/include/VUSolid.hh"

class Geometry {
public:
  void AddNode(const VUSolid* s) { solids_.push_back(s); }
  VUSolid::EnumInside Inside (const UVector3 &aPoint) const;
  const VUSolid* InsideWhat(const UVector3 &aPoint) const;
  double  SafetyFromInside ( const UVector3 &aPoint, bool aAccurate=false) const;
  double  SafetyFromOutside( const UVector3 &aPoint, bool aAccurate=false) const;
private:
  std::vector<const VUSolid*> solids_;
};
#endif
