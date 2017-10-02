#ifndef WITH_USOLIDS
#include "SimpleGeom.h"

VUSolid::EnumInside VUSolid::Inside (const UVector3 &aPoint) const
{
  double z = aPoint.Z();
  if (zmin_ < z && z < zmax_)
  {
    double r = aPoint.Perp();
    if (rin_ < r && r < rout_) return eInside;
  }
  return eOutside;
}

double VUSolid::SafetyFromInside (const UVector3 &aPoint, bool) const
{
  double r = aPoint.Perp();
  if (rin_ < r && r < rout_)
  {
    double z = aPoint.Z();
    if (zmin_ < z && z < zmax_)
    {
      double dr = std::min(r - rin_, rout_ - r);
      return std::max(dr, std::min(z - zmin_, zmax_ - z));
    }
  }
  return 0.0;
}

double VUSolid::SafetyFromOutside(const UVector3 &aPoint, bool) const
{
  double r = aPoint.Perp(), z = aPoint.Z(), dr, dz;
  if      (r < rin_)  dr = rin_ - r;
  else if (r > rout_) dr = r - rout_;
  else                dr = 0;
  if      (z < zmin_) dz = zmin_ - z;
  else if (z > zmax_) dz = z - zmax_;
  else                dz = 0;
  // XXXX return std::hypot(dr, dz);
  return std::max(dr, dz);
}

double VUSolid::SafetyFromOutsideDr(const UVector3 &aPoint, double ooaCtgTheta, bool) const
{
  double r = aPoint.Perp(), z = aPoint.Z(), dr, dz;
  if      (r < rin_)  dr = rin_ - r;
  else if (r > rout_) dr = r - rout_;
  else                dr = 0;
  if      (z < zmin_) dz = zmin_ - z;
  else if (z > zmax_) dz = z - zmax_;
  else                dz = 0;

  return std::max(dr, dz * ooaCtgTheta);
}

bool VUSolid::Normal(const UVector3& aPoint, UVector3 &aNormal) const
{
  const double tol = 1e-6;
  const double rho = aPoint.Perp();
  const double zoo = aPoint.At(2);

  bool nz_is_min = std::abs(zoo - zmin_) < std::abs(zoo - zmax_);
  double dz = nz_is_min ? std::abs(zoo - zmin_) : std::abs(zoo - zmax_);

  bool nr_is_in  = std::abs(rho - rin_) < std::abs(rho - rout_);
  double dr = nr_is_in  ? std::abs(rho - rin_) : std::abs(rho - rout_);

  if (dr < tol && dz < tol)
  {
     aNormal.At(0) = aPoint.At(0) / rho * (nr_is_in ? - Config::OOSqrt2 : Config::OOSqrt2);
     aNormal.At(1) = aPoint.At(1) / rho * (nr_is_in ? - Config::OOSqrt2 : Config::OOSqrt2);
     aNormal.At(2) = nz_is_min ? - Config::OOSqrt2 : Config::OOSqrt2;
  }
  else if (dr < dz)
  {
    if (nr_is_in)
      aNormal = UVector3(-aPoint.At(0)/rho, -aPoint.At(1)/rho, 0.0);
    else
      aNormal = UVector3( aPoint.At(0)/rho,  aPoint.At(1)/rho, 0.0);
  }
  else
  {
    aNormal = UVector3(0.0, 0.0, nz_is_min ? -1.0 : 1.0);
  }

  return dz < tol || dr < tol;
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
