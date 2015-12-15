#ifndef _propagation_mplex_
#define _propagation_mplex_

#include "Track.h"
#include "Matrix.h"

void propagateLineToRMPlex(const MPlexLS &psErr,  const MPlexLV& psPar,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                                 MPlexLS &outErr,       MPlexLV& outPar);

void propagateHelixToRMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const MPlexHV& msPar,
			          MPlexLS &outErr,       MPlexLV& outPar);

void propagateHelixToRMPlex(const MPlexLS& inErr,  const MPlexLV& inPar,
                            const MPlexQI& inChg,  const float    r,
			          MPlexLS& outErr,       MPlexLV& outPar,
                            const int      N_proc);

//inline?
inline void computeJacobianSimple(int n, MPlexLL& errorProp, 
				  float s, float k, float p, float pxin, float pyin, float pzin, 
				  float TP, float cosTP, float sinTP);

void helixAtRFromIterative(const MPlexLV& inPar, const MPlexQI& inChg, 
			         MPlexLV& outPar, const MPlexQF &msRad, 
			         MPlexLL& errorProp, bool useSimpleJac);

void helixAtRFromIntersection(const MPlexLV& inPar, const MPlexQI& inChg, 
                                    MPlexLV& outPar, const MPlexQF &msRad, 
   			            MPlexLL& errorProp);

void applyMaterialEffects(const MPlexQF &hitsRl, const MPlexQF& hitsXi, 
			  MPlexLS &outErr, MPlexLV& outPar);

inline float getRlVal(const float r, const float zin) {
  float z = fabs(zin);
  //pixel
  if (r<7) return 0.035;
  if (r<9) {
    if (z<7) return 0.020;
    else     return 0.020;
  }
  if (r<20) {
    if (z<12) return 0.020;
    else      return 0.025;
  }
  //TIB
  if (r<30) {
    if (z<22) return 0.025;
    else      return 0.035;
  }
  if (r<38) {
    if (z<27) return 0.025;
    else      return 0.040;
  }
  if (r<46) {
    if (z<45) return 0.040;
    else      return 0.060;
  }
  if (r<55) {
    if (z<50) return 0.040;
    else      return 0.050;
  }
  //TOB
  if (r<65) {
    if (z<17)      return 0.020;
    else if (z<70) return 0.040;
    else           return 0.060;
  }
  if (r<75) {
    if (z<17)      return 0.010;
    else if (z<70) return 0.025;
    else           return 0.030;
  }
  if (r<82) {
    if (z<17)      return 0.015;
    else if (z<70) return 0.040;
    else           return 0.045;
  }
  if (r<90) {
    if (z<17)      return 0.015;
    else if (z<70) return 0.040;
    else           return 0.040;
  }
  if (r<100) {
    if (z<17)      return 0.015;
    else if (z<70) return 0.040;
    else           return 0.040;
  }
  if (r<102) {
    if (z<17)      return 0.015;
    else if (z<70) return 0.040;
    else           return 0.040;
  }
  return 0.;
}

inline float getXiVal(const float r, const float zin) {
  float z = fabs(zin);
  //pixel
  if (r<7) return 0.08e-03;
  if (r<9) {
    if (z<7) return 0.03e-03;
    else     return 0.04e-03;
  }
  if (r<20) {
    if (z<12) return 0.03e-03;
    else      return 0.04e-03;
  }
  //TIB
  if (r<30) {
    if (z<22) return 0.06e-03;
    else      return 0.08e-03;
  }
  if (r<38) {
    if (z<27) return 0.04e-03;
    else      return 0.06e-03;
  }
  if (r<46) {
    if (z<45) return 0.07e-03;
    else      return 0.12e-03;
  }
  if (r<55) {
    if (z<50) return 0.07e-03;
    else      return 0.11e-03;
  }
  //TOB
  if (r<65) {
    if (z<17)      return 0.05e-03;
    else if (z<70) return 0.07e-03;
    else           return 0.12e-03;
  }
  if (r<75) {
    if (z<17)      return 0.02e-03;
    else if (z<70) return 0.04e-03;
    else           return 0.05e-03;
  }
  if (r<82) {
    if (z<17)      return 0.02e-03;
    else if (z<70) return 0.07e-03;
    else           return 0.07e-03;
  }
  if (r<90) {
    if (z<17)      return 0.02e-03;
    else if (z<70) return 0.07e-03;
    else           return 0.07e-03;
  }
  if (r<100) {
    if (z<17)      return 0.02e-03;
    else if (z<70) return 0.07e-03;
    else           return 0.07e-03;
  }
  if (r<102) {
    if (z<17)      return 0.02e-03;
    else if (z<70) return 0.07e-03;
    else           return 0.07e-03;
  }
  return 0.;
}

#endif
