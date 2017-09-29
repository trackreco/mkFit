#ifndef _simplegeom_
#define _simplegeom_
#include <vector>
#include "Matrix.h"

// MT XXXX: Should SVector3 base class really be the double
//          version (probably is float right now).

#ifndef WITH_USOLIDS
class UVector3 : public SVector3
{
  public:
  UVector3(float x, float y, float z) : SVector3(x, y, z) {}
  UVector3() : SVector(0.0f,0.0f,0.0f) {}

  inline double X() const { return At(0); }
  inline double Y() const { return At(1); }
  inline double Z() const { return At(2); }

  inline double Dot(const UVector3&) const;
  inline double Mag2() const;
  inline double Mag() const;
  inline double Perp2() const;
  inline double Perp() const;
  double Normalize();
};

inline double UVector3::Dot(const UVector3& p) const
{
  return At(0) * p.At(0) + At(1) * p.At(1) + At(2) * p.At(2);
}

inline double UVector3::Mag2() const
{
  return At(0) * At(0) + At(1) * At(1) + At(2) * At(2);
}

inline double UVector3::Mag() const
{
  return std::sqrt(Mag2());
}

inline double UVector3::Perp2() const
{
  return At(0) * At(0) + At(1) * At(1);
}

inline double UVector3::Perp() const
{
  return std::sqrt(Perp2());
}

class VUSolid
{
public:
  enum EnumInside { eInside=0, eSurface=1, eOutside=2 }; 

  VUSolid(double rin, double rout, double zmin, double zmax,
          bool is_barrel, bool is_outer) :
    rin_(rin), rout_(rout), zmin_(zmin), zmax_(zmax),
    is_barrel_(is_barrel), is_outer_(is_outer)
  {}

  EnumInside Inside (const UVector3 &aPoint) const;
  double SafetyFromInside (const UVector3 &aPoint, bool aAccurate=false) const;
  double SafetyFromOutside(const UVector3 &aPoint, bool aAccurate=false) const;
  double SafetyFromOutsideDr(const UVector3 &aPoint, double ooaCtgTheta, bool aAccurate=false) const;
  bool   Normal(const UVector3& aPoint, UVector3 &aNormal) const;

  VUSolid* Clone() const {return new VUSolid(rin_, rout_, zmin_, zmax_, is_barrel_, is_outer_);}

  double rin_,  rout_;
  double zmin_, zmax_;
  bool   is_barrel_;
  bool   is_outer_;
};

#endif
#endif
