#ifndef MkBase_h
#define MkBase_h

#include "Matrix.h"

#include "PropagationMPlex.h"

//==============================================================================
// MkBase
//==============================================================================

struct MkBase
{
  MPlexLS    Err[2];
  MPlexLV    Par[2];

  MPlexQI    Chg;

  static constexpr int iC = 0; // current
  static constexpr int iP = 1; // propagated

  float getPar(int itrack, int i, int par) const { return Par[i].ConstAt(itrack, 0, par); }

  //----------------------------------------------------------------------------

  MkBase() {}

  //----------------------------------------------------------------------------

  void PropagateTracksToR(float r, const int N_proc, const PropagationFlags pf)
  {
    MPlexQF msRad;
#pragma simd
    for (int n = 0; n < NN; ++n)
    {
      msRad.At(n, 0, 0) = r;
    }

    propagateHelixToRMPlex(Err[iC], Par[iC], Chg, msRad,
                           Err[iP], Par[iP], N_proc, pf);
  }

  void PropagateTracksToHitR(const MPlexHV& par, const int N_proc, const PropagationFlags pf)
  {
    MPlexQF msRad;
#pragma simd
    for (int n = 0; n < NN; ++n)
    {
      msRad.At(n, 0, 0) = std::hypot(par.ConstAt(n, 0, 0), par.ConstAt(n, 1, 0));
    }

    propagateHelixToRMPlex(Err[iC], Par[iC], Chg, msRad,
                           Err[iP], Par[iP], N_proc, pf);
  }

  //----------------------------------------------------------------------------

  void PropagateTracksToZ(float z, const int N_proc, const PropagationFlags pf)
  {
    MPlexQF msZ;
#pragma simd
    for (int n = 0; n < NN; ++n)
    {
      msZ.At(n, 0, 0) = z;
    }

    propagateHelixToZMPlex(Err[iC], Par[iC], Chg, msZ,
                           Err[iP], Par[iP], N_proc, pf);
  }

  void PropagateTracksToHitZ(const MPlexHV& par, const int N_proc, const PropagationFlags pf)
  {
    MPlexQF msZ;
#pragma simd
    for (int n = 0; n < NN; ++n)
    {
      msZ.At(n, 0, 0) = par.ConstAt(n, 2, 0);
    }

    propagateHelixToZMPlex(Err[iC], Par[iC], Chg, msZ,
                           Err[iP], Par[iP], N_proc, pf);
  }

};

#endif
