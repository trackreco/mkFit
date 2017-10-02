#ifndef _geometry_
#define _geometry_

#include <vector>
#include "Matrix.h"

#ifdef WITH_USOLIDS
#include "USolids/include/VUSolid.hh"
#else
#include "SimpleGeom.h"
#endif

class TrackerInfo;

// XXXX MT: What it the purpose of zPlanes?
// I guess this is half-baked and we actually need two add layer
// methods, like AddBarellLayer(), AddEndcapLayer().
// And probably something even more descriptive, so we can also
// build navigation LUTs (or something).

class Geometry
{
public:
  ~Geometry() { for (auto& s : solids_) { delete s; }}
  void AddLayer(const VUSolid* s) { solids_.push_back(s); }
  void AddLayer(const VUSolid* s, const float radius, const float zPlane) { 
    solids_.push_back(s); radii_.push_back(radius); zPlanes_.push_back(zPlane);
  }
  int CountLayers() const { return solids_.size(); }
  const VUSolid* Layer(int l) const { return solids_.at(l); }
  int   LayerOfSolid(const VUSolid *s) const;
  float Radius(int l) const { return radii_.at(l); }
  float zPlane(int l) const { return zPlanes_.at(l); }
  VUSolid::EnumInside Inside (const UVector3 &aPoint) const;
  int LayerIndex(const UVector3 &aPoint) const;
  const VUSolid* InsideWhat(const UVector3 &aPoint) const;
  double  SafetyFromInside   ( const UVector3 &aPoint, bool aAccurate=false) const;
  double  SafetyFromOutside  ( const UVector3 &aPoint, bool aAccurate=false) const;
  double  SafetyFromOutsideDr( const UVector3 &aPoint, double ooaCtgTheta, int skip_layer, int &layer, bool aAccurate=false) const;
  Geometry clone() const;

  void BuildFromTrackerInfo(const TrackerInfo& tracker_info);

private:
  std::vector<const VUSolid*> solids_;
  std::vector<float> radii_;
  std::vector<float> zPlanes_;
};

#endif
