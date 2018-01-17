#ifndef _propagation_mplex_
#define _propagation_mplex_

#include "Matrix.h"

inline void squashPhiMPlex(MPlexLV& par, const int N_proc)
{
  #pragma simd
  for (int n = 0; n < NN; ++n) {
    if (par(n, 4, 0) >= Config::PI) par(n, 4, 0) -= Config::TwoPI;
    if (par(n, 4, 0) < -Config::PI) par(n, 4, 0) += Config::TwoPI;
  }
}

inline void squashPhiMPlexGeneral(MPlexLV& par, const int N_proc)
{
  #pragma simd
  for (int n = 0; n < NN; ++n) {
    par(n, 4, 0) -= std::floor(0.5f*Config::InvPI*(par(n, 4, 0)+Config::PI)) * Config::TwoPI;
  }
}

void propagateLineToRMPlex(const MPlexLS &psErr,  const MPlexLV& psPar,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                                 MPlexLS &outErr,       MPlexLV& outPar,
                           const int      N_proc);

void propagateHelixToRMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const MPlexQF& msRad,
                                  MPlexLS &outErr,       MPlexLV& outPar,
                            const int      N_proc, const PropagationFlags pflags);

void helixAtRFromIterativeCCSFullJac(const MPlexLV& inPar, const MPlexQI& inChg, const MPlexQF &msRad,
                                           MPlexLV& outPar,      MPlexLL& errorProp,
                                     const int      N_proc);

void helixAtRFromIterativeCCS(const MPlexLV& inPar,  const MPlexQI& inChg, const MPlexQF &msRad,
                                    MPlexLV& outPar,       MPlexLL& errorProp,
                              const int      N_proc, const PropagationFlags pflags);

void propagateHelixToZMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const MPlexQF& msZ,
                                  MPlexLS &outErr,       MPlexLV& outPar,
                            const int      N_proc, const PropagationFlags pflags);

void helixAtZ(const MPlexLV& inPar,  const MPlexQI& inChg, const MPlexQF &msZ,
                    MPlexLV& outPar,       MPlexLL& errorProp,
              const int      N_proc, const PropagationFlags pflags);

void applyMaterialEffects(const MPlexQF &hitsRl, const MPlexQF& hitsXi, 
                                MPlexLS &outErr,       MPlexLV& outPar,
                          const int      N_proc);

#endif
