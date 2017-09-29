#ifndef MkBase_h
#define MkBase_h

#include "Matrix.h"

#include "PropagationMPlex.h"

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

  void PropagateTracksToR(float R, const int N_proc)
  {
    propagateHelixToRMPlex(Err[iC], Par[iC], Chg, R,
                           Err[iP], Par[iP], N_proc);
  }

  void PropagateTracksToZ(float Z, const int N_proc)
  {
    propagateHelixToZMPlex(Err[iC], Par[iC], Chg, Z,
                           Err[iP], Par[iP], N_proc);
  }

};

#endif
