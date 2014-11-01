#ifndef WITH_USOLIDS
#include "SimpleGeom.h"

VUSolid::EnumInside VUSolid::Inside (const UVector3 &aPoint) const
{
  double r = aPoint.Perp();
  if (rin_ < r && r < rout_) return eInside;
  return eOutside;
}

double VUSolid::SafetyFromInside (const UVector3 &aPoint, bool) const
{
  double r = aPoint.Perp();
  if (rin_ < r && r < rout_) return std::min(r - rin_, rout_ - r);
  return 0.0;
}

double VUSolid::SafetyFromOutside(const UVector3 &aPoint, bool) const
{
  double r = aPoint.Perp();
  if (rin_ < r && r < rout_) return 0.0;
  return std::min(std::abs(r - rin_), std::abs(rout_ - r));
}

bool VUSolid::Normal( const UVector3& aPoint, UVector3 &aNormal ) const
{
  double rho = aPoint.Perp();
  aNormal = UVector3(aPoint.At(0)/rho, aPoint.At(1)/rho, 0.0);
  return true;
}

double UVector3::Normalize()
{
  // Normalize to unit. Return normalization factor.
  double  mag = Mag2();
  if (mag == 0.0) return mag;;
  mag = std::sqrt(mag);
  At(0) /= mag;
  At(1) /= mag;
  At(2) /= mag;
  return mag;
}
#endif
