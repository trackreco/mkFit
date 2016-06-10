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
                            const int      N_proc);

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
                              const int      N_proc);

void helixAtRFromIterative(const MPlexLV& inPar, const MPlexQI& inChg, 
			         MPlexLV& outPar, const MPlexQF &msRad, 
			         MPlexLL& errorProp, bool useSimpleJac,
                           const int      N_proc);

void applyMaterialEffects(const MPlexQF &hitsRl, const MPlexQF& hitsXi, 
                                MPlexLS &outErr, MPlexLV& outPar,
                          const int      N_proc);

inline float getRlVal(const float r, const float zin) {
  float z = std::abs(zin);
  //pixel
  if (r<6) return 0.026f;
  if (r<9) {
    if (z<7) return 0.020f;
    else     return 0.020f;
  }
  if (r<20) {
    if (z<12) return 0.020f;
    else      return 0.025f;
  }
  //TIB
  if (r<30) {
    if (z<22) return 0.025f;
    else      return 0.035f;
  }
  if (r<38) {
    if (z<27) return 0.023f;
    else      return 0.040f;
  }
  if (r<46) {
    if (z<45) return 0.040f;
    else      return 0.060f;
  }
  if (r<55) {
    if (z<50) return 0.036f;
    else      return 0.050f;
  }
  //TOB
  if (r<65) {
    if (z<17)      return 0.022f;
    else if (z<70) return 0.040f;
    else           return 0.060f;
  }
  if (r<75) {
    if (z<17)      return 0.012f;
    else if (z<70) return 0.025f;
    else           return 0.030f;
  }
  if (r<82) {
    if (z<17)      return 0.015f;
    else if (z<70) return 0.040f;
    else           return 0.045f;
  }
  if (r<90) {
    if (z<17)      return 0.015f;
    else if (z<70) return 0.040f;
    else           return 0.040f;
  }
  if (r<100) {
    if (z<17)      return 0.015f;
    else if (z<70) return 0.040f;
    else           return 0.040f;
  }
  if (r<102) {
    if (z<17)      return 0.015f;
    else if (z<70) return 0.040f;
    else           return 0.040f;
  }
  return 0.0f;
}

inline float getXiVal(const float r, const float zin) {
  float z = std::abs(zin);
  //pixel
  if (r<6) return 0.054e-03f;
  if (r<9) {
    if (z<7) return 0.037e-03f;
    else     return 0.04e-03f;
  }
  if (r<20) {
    if (z<12) return 0.037e-03f;
    else      return 0.04e-03f;
  }
  //TIB
  if (r<30) {
    if (z<22) return 0.057e-03f;
    else      return 0.08e-03f;
  }
  if (r<38) {
    if (z<27) return 0.048e-03f;
    else      return 0.06e-03f;
  }
  if (r<46) {
    if (z<45) return 0.085e-03f;
    else      return 0.12e-03f;
  }
  if (r<55) {
    if (z<50) return 0.075e-03f;
    else      return 0.11e-03f;
  }
  //TOB
  if (r<65) {
    if (z<17)      return 0.056e-03f;
    else if (z<70) return 0.07e-03f;
    else           return 0.12e-03f;
  }
  if (r<75) {
    if (z<17)      return 0.024e-03f;
    else if (z<70) return 0.04e-03f;
    else           return 0.05e-03f;
  }
  if (r<82) {
    if (z<17)      return 0.034e-03f;
    else if (z<70) return 0.07e-03f;
    else           return 0.07e-03f;
  }
  if (r<90) {
    if (z<17)      return 0.033e-03f;
    else if (z<70) return 0.07e-03f;
    else           return 0.07e-03f;
  }
  if (r<100) {
    if (z<17)      return 0.033e-03f;
    else if (z<70) return 0.07e-03f;
    else           return 0.07e-03f;
  }
  if (r<102) {
    if (z<17)      return 0.033e-03f;
    else if (z<70) return 0.07e-03f;
    else           return 0.07e-03f;
  }
  return 0.0f;
}

#endif
