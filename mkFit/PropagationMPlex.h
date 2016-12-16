#ifndef _propagation_mplex_
#define _propagation_mplex_

#include "Track.h"
#include "Matrix.h"

void propagateLineToRMPlex(const MPlexLS &psErr,  const MPlexLV& psPar,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                                 MPlexLS &outErr,       MPlexLV& outPar,
                           const int      N_proc);

void propagateHelixToRMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const MPlexHV& msPar,
                                  MPlexLS &outErr,       MPlexLV& outPar,
                            const int      N_proc, const bool useParamBfield = false);

void propagateHelixToRMPlex(const MPlexLS& inErr,  const MPlexLV& inPar,
                            const MPlexQI& inChg,  const float    r,
                                  MPlexLS& outErr,       MPlexLV& outPar,
                            const int      N_proc);

void helixAtRFromIterativeCCSFullJac(const MPlexLV& inPar, const MPlexQI& inChg,
                                           MPlexLV& outPar, const MPlexQF &msRad,
                                           MPlexLL& errorProp,
                                     const int      N_proc);

void helixAtRFromIterativeCCS(const MPlexLV& inPar, const MPlexQI& inChg,
                                    MPlexLV& outPar, const MPlexQF &msRad,
                                    MPlexLL& errorProp,
                              const int      N_proc, const bool useParamBfield = false);

void helixAtRFromIterative(const MPlexLV& inPar, const MPlexQI& inChg, 
                                 MPlexLV& outPar, const MPlexQF &msRad, 
                                 MPlexLL& errorProp,
                           const int      N_proc, const bool useParamBfield = false);

void propagateHelixToZMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const MPlexHV& msPar,
                                  MPlexLS &outErr,       MPlexLV& outPar,
                            const int      N_proc, const bool useParamBfield = false);

void propagateHelixToZMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const float    z,
                                  MPlexLS &outErr,       MPlexLV& outPar,
                            const int      N_proc);

void helixAtZ(const MPlexLV& inPar,  const MPlexQI& inChg,
                    MPlexLV& outPar, const MPlexQF &msZ,
                    MPlexLL& errorProp,
              const int      N_proc, const bool useParamBfield = false);

void applyMaterialEffects(const MPlexQF &hitsRl, const MPlexQF& hitsXi, 
                                MPlexLS &outErr, MPlexLV& outPar,
                          const int      N_proc);

inline float getRlVal(const float r, const float zin) {
  float z = std::abs(zin);
  //pixel barrel
  if (r<6 && z<30) return 0.026f;
  if (r<9 && z<30) {
    if (z<7) return 0.020f;
    else     return 0.020f;
  }
  if (r<20 && z<30) {
    if (z<12) return 0.020f;
    else      return 0.025f;
  }
  //pixel endcap
  if (r<20 && z>30) {
    if (z<40) return 0.109f;
    else      return 0.068f;
  }
  //TIB
  if (r<30 && z<70) {
    if (z<22) return 0.025f;
    else      return 0.035f;
  }
  if (r<38 && z<70) {
    if (z<27) return 0.023f;
    else      return 0.040f;
  }
  if (r<46 && z<70) {
    if (z<45) return 0.040f;
    else      return 0.060f;
  }
  if (r<55 && z<70) {
    if (z<50) return 0.036f;
    else      return 0.050f;
  }
  //TID
  if (r<55 && z>70 && z<120) {
    if (z<81)       return 0.180f;
    else if (z<94)  return 0.090f;
    else if (z<99)  return 0.070f;
    else if (z<103) return 0.096f;
    else            return 0.048f;
  }
  //TOB
  if (r<65 && z<120) {
    if (z<17)      return 0.022f;
    else if (z<70) return 0.040f;
    else           return 0.060f;
  }
  if (r<75 && z<120) {
    if (z<17)      return 0.012f;
    else if (z<70) return 0.025f;
    else           return 0.030f;
  }
  if (r<82 && z<120) {
    if (z<17)      return 0.015f;
    else if (z<70) return 0.040f;
    else           return 0.045f;
  }
  if (r<90 && z<120) {
    if (z<17)      return 0.015f;
    else if (z<70) return 0.040f;
    else           return 0.040f;
  }
  if (r<100 && z<120) {
    if (z<17)      return 0.015f;
    else if (z<70) return 0.040f;
    else           return 0.040f;
  }
  if (r<102 && z<120) {
    if (z<17)      return 0.015f;
    else if (z<70) return 0.040f;
    else           return 0.040f;
  }
  //TEC
  if (z>120) {
    if (z<128) {
      if (r<55)      return 0.103f;
      else if (r<80) return 0.145f;
      else           return 0.279f;
    }
    if (z<132) {
      if (r<45)      return 0.103f;
      else if (r<70) return 0.145f;
      else           return 0.279f;
    }
    if (z<136) {
      if (r<55)      return 0.103f;
      else if (r<80) return 0.145f;
      else           return 0.279f;
    }
    if (z<138) {
      if (r<45)      return 0.103f;
      else if (r<70) return 0.145f;
      else           return 0.279f;
    }

    if (z<143) {
      if (r<35)      return 0.036f;
      else if (r<55) return 0.071f;
      else if (r<80) return 0.036f;
      else           return 0.071f;
    }
    if (z<146) {
      if (r<45)      return 0.036f;
      else           return 0.071f;
    }

    if (z<150) {
      if (r<35)      return 0.036f;
      else if (r<55) return 0.071f;
      else if (r<80) return 0.036f;
      else           return 0.071f;
    }
    if (z<153) {
      if (r<45)      return 0.036f;
      else           return 0.071f;
    }

    if (z<157) {
      if (r<35)      return 0.027f;
      else if (r<55) return 0.054f;
      else if (r<80) return 0.027f;
      else           return 0.054f;
    }
    if (z<160) {
      if (r<45)      return 0.027f;
      else           return 0.054f;
    }

    if (z<164) {
      if (r<35)      return 0.027f;
      else if (r<55) return 0.054f;
      else if (r<80) return 0.027f;
      else           return 0.054f;
    }
    if (z<167) {
      if (r<45)      return 0.027f;
      else           return 0.054f;
    }

    if (z<170) {
      if (r<55)      return 0.047f;
      else if (r<80) return 0.024f;
      else           return 0.047f;
    }
    if (z<175) {
      if (r<45)      return 0.024f;
      else           return 0.047f;
    }

    if (z<178) {
      if (r<55)      return 0.047f;
      else if (r<80) return 0.024f;
      else           return 0.047f;
    }
    if (z<180) {
      if (r<45)      return 0.024f;
      else           return 0.047f;
    }

    if (z<185) {
      if (r<55)      return 0.051f;
      else if (r<80) return 0.025f;
      else           return 0.051f;
    }
    if (z<188) {
      if (r<45)      return 0.025f;
      else           return 0.051f;
    }

    if (z<192) {
      if (r<55)      return 0.051f;
      else if (r<80) return 0.025f;
      else           return 0.051f;
    }
    if (z<195) {
      if (r<45)      return 0.025f;
      else           return 0.051f;
    }

    if (z<202) {
      if (r<55)      return 0.053f;
      else if (r<80) return 0.026f;
      else           return 0.053f;
    }
    if (z<205) {
      if (r<45)      return 0.026f;
      else           return 0.053f;
    }

    if (z<209) {
      if (r<55)      return 0.053f;
      else if (r<80) return 0.026f;
      else           return 0.053f;
    }
    if (z<212) {
      if (r<45)      return 0.026f;
      else           return 0.053f;
    }

    if (z<221) {
      if (r<55)      return 0.047f;
      else if (r<80) return 0.024f;
      else           return 0.047f;
    }
    if (z<224)       return 0.047f;

    if (z<228) {
      if (r<55)      return 0.047f;
      else if (r<80) return 0.024f;
      else           return 0.047f;
    }
    if (z<232)       return 0.047f;

    if (z<242) {
      if (r<55)      return 0.045f;
      else if (r<80) return 0.023f;
      else           return 0.045f;
    }
    if (z<244)       return 0.045f;

    if (z<249) {
      if (r<55)      return 0.045f;
      else if (r<80) return 0.023f;
      else           return 0.045f;
    }
    if (z<252)       return 0.045f;

    if (z<265) {
      if (r<80)      return 0.045f;
      else           return 0.023f;
    }
    if (z<252)       return 0.045f;

    if (z<270) {
      if (r<80)      return 0.045f;
      else           return 0.023f;
    }
    if (z<280)       return 0.045f;
  }
  return 0.0f;
}

inline float getXiVal(const float r, const float zin) {
  float z = std::abs(zin);
  //pixel barrel
  if (r<6 && z<30) return 0.054e-03f;
  if (r<9 && z<30) {
    if (z<7) return 0.037e-03f;
    else     return 0.04e-03f;
  }
  if (r<20 && z<30) {
    if (z<12) return 0.037e-03f;
    else      return 0.04e-03f;
  }
  //pixel endcap
  if (r<20 && z>30) {
    if (z<40) return 0.20e-03f;
    else      return 0.13e-03f;
  }
  //TIB
  if (r<30 && z<70) {
    if (z<22) return 0.057e-03f;
    else      return 0.08e-03f;
  }
  if (r<38 && z<70) {
    if (z<27) return 0.048e-03f;
    else      return 0.06e-03f;
  }
  if (r<46 && z<70) {
    if (z<45) return 0.085e-03f;
    else      return 0.12e-03f;
  }
  if (r<55 && z<70) {
    if (z<50) return 0.075e-03f;
    else      return 0.11e-03f;
  }
  //TID
  if (r<55 && z>70 && z<120) {
    if (z<81)       return 0.34e-03f;
    else if (z<94)  return 0.17e-03f;
    else if (z<99)  return 0.07e-03f;
    else if (z<103) return 0.22e-03f;
    else            return 0.11e-03f;
  }
  //TOB
  if (r<65 && z<120) {
    if (z<17)      return 0.056e-03f;
    else if (z<70) return 0.07e-03f;
    else           return 0.12e-03f;
  }
  if (r<75 && z<120) {
    if (z<17)      return 0.024e-03f;
    else if (z<70) return 0.04e-03f;
    else           return 0.05e-03f;
  }
  if (r<82 && z<120) {
    if (z<17)      return 0.034e-03f;
    else if (z<70) return 0.07e-03f;
    else           return 0.07e-03f;
  }
  if (r<90 && z<120) {
    if (z<17)      return 0.033e-03f;
    else if (z<70) return 0.07e-03f;
    else           return 0.07e-03f;
  }
  if (r<100 && z<120) {
    if (z<17)      return 0.033e-03f;
    else if (z<70) return 0.07e-03f;
    else           return 0.07e-03f;
  }
  if (r<102 && z<120) {
    if (z<17)      return 0.033e-03f;
    else if (z<70) return 0.07e-03f;
    else           return 0.07e-03f;
  }
  //TEC
  if (z>120) {
    if (z<128) {
      if (r<80)      return 0.27e-03f;
      else           return 0.48e-03f;
    }
    if (z<132) {
      if (r<45)      return 0.24e-03f;
      else if (r<70) return 0.39e-03f;
      else           return 0.55e-03f;
    }
    if (z<136) {
      if (r<80)      return 0.27e-03f;
      else           return 0.48e-03f;
    }
    if (z<138) {
      if (r<45)      return 0.24e-03f;
      else if (r<70) return 0.39e-03f;
      else           return 0.55e-03f;
    }

    if (z<143) {
      if (r<35)      return 0.08e-03f;
      else if (r<55) return 0.16e-03f;
      else if (r<80) return 0.08e-03f;
      else           return 0.16e-03f;
    }
    if (z<146) {
      if (r<45)      return 0.08e-03f;
      else           return 0.16e-03f;
    }

    if (z<150) {
      if (r<35)      return 0.08e-03f;
      else if (r<55) return 0.16e-03f;
      else if (r<80) return 0.08e-03f;
      else           return 0.16e-03f;
    }
    if (z<153) {
      if (r<45)      return 0.08e-03f;
      else           return 0.16e-03f;
    }

    if (z<157) {
      if (r<35)      return 0.06e-03f;
      else if (r<55) return 0.12e-03f;
      else if (r<80) return 0.06e-03f;
      else           return 0.12e-03f;
    }
    if (z<160) {
      if (r<45)      return 0.06e-03f;
      else           return 0.12e-03f;
    }

    if (z<164) {
      if (r<35)      return 0.06e-03f;
      else if (r<55) return 0.12e-03f;
      else if (r<80) return 0.06e-03f;
      else           return 0.12e-03f;
    }
    if (z<167) {
      if (r<45)      return 0.06e-03f;
      else           return 0.12e-03f;
    }

    if (z<170) {
      if (r<55)      return 0.11e-03f;
      else if (r<80) return 0.05e-03f;
      else           return 0.11e-03f;
    }
    if (z<175) {
      if (r<45)      return 0.05e-03f;
      else           return 0.11e-03f;
    }

    if (z<178) {
      if (r<55)      return 0.11e-03f;
      else if (r<80) return 0.05e-03f;
      else           return 0.11e-03f;
    }
    if (z<180) {
      if (r<45)      return 0.05e-03f;
      else           return 0.11e-03f;
    }

    if (z<185) {
      if (r<55)      return 0.12e-03f;
      else if (r<80) return 0.06e-03f;
      else           return 0.12e-03f;
    }
    if (z<188) {
      if (r<45)      return 0.06e-03f;
      else           return 0.12e-03f;
    }

    if (z<192) {
      if (r<55)      return 0.12e-03f;
      else if (r<80) return 0.06e-03f;
      else           return 0.12e-03f;
    }
    if (z<195) {
      if (r<45)      return 0.06e-03f;
      else           return 0.12e-03f;
    }

    if (z<202) {
      if (r<55)      return 0.12e-03f;
      else if (r<80) return 0.06e-03f;
      else           return 0.12e-03f;
    }
    if (z<205) {
      if (r<45)      return 0.06e-03f;
      else           return 0.12e-03f;
    }

    if (z<209) {
      if (r<55)      return 0.12e-03f;
      else if (r<80) return 0.06e-03f;
      else           return 0.12e-03f;
    }
    if (z<212) {
      if (r<45)      return 0.06e-03f;
      else           return 0.12e-03f;
    }

    if (z<221) {
      if (r<55)      return 0.11e-03f;
      else if (r<80) return 0.05e-03f;
      else           return 0.11e-03f;
    }
    if (z<224)       return 0.11e-03f;

    if (z<228) {
      if (r<55)      return 0.11e-03f;
      else if (r<80) return 0.05e-03f;
      else           return 0.11e-03f;
    }
    if (z<232)       return 0.11e-03f;

    if (z<242) {
      if (r<55)      return 0.10e-03f;
      else if (r<80) return 0.05e-03f;
      else           return 0.10e-03f;
    }
    if (z<244)       return 0.10e-03f;

    if (z<249) {
      if (r<55)      return 0.10e-03f;
      else if (r<80) return 0.05e-03f;
      else           return 0.10e-03f;
    }
    if (z<252)       return 0.10e-03f;

    if (z<265) {
      if (r<80)      return 0.05e-03f;
      else           return 0.10e-03f;
    }
    if (z<252)       return 0.10e-03f;

    if (z<270) {
      if (r<80)      return 0.05e-03f;
      else           return 0.10e-03f;
    }
    if (z<280)       return 0.10e-03f;
  }
  return 0.0f;
}

#endif
