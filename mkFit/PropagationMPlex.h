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
				  float k, float TP, float cosTP, float sinTP);

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
  if (r<6) return 0.026;
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
    if (z<27) return 0.023;
    else      return 0.040;
  }
  if (r<46) {
    if (z<45) return 0.040;
    else      return 0.060;
  }
  if (r<55) {
    if (z<50) return 0.036;
    else      return 0.050;
  }
  //TOB
  if (r<65) {
    if (z<17)      return 0.022;
    else if (z<70) return 0.040;
    else           return 0.060;
  }
  if (r<75) {
    if (z<17)      return 0.012;
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
  if (r<6) return 0.054e-03;
  if (r<9) {
    if (z<7) return 0.037e-03;
    else     return 0.04e-03;
  }
  if (r<20) {
    if (z<12) return 0.037e-03;
    else      return 0.04e-03;
  }
  //TIB
  if (r<30) {
    if (z<22) return 0.057e-03;
    else      return 0.08e-03;
  }
  if (r<38) {
    if (z<27) return 0.048e-03;
    else      return 0.06e-03;
  }
  if (r<46) {
    if (z<45) return 0.085e-03;
    else      return 0.12e-03;
  }
  if (r<55) {
    if (z<50) return 0.075e-03;
    else      return 0.11e-03;
  }
  //TOB
  if (r<65) {
    if (z<17)      return 0.056e-03;
    else if (z<70) return 0.07e-03;
    else           return 0.12e-03;
  }
  if (r<75) {
    if (z<17)      return 0.024e-03;
    else if (z<70) return 0.04e-03;
    else           return 0.05e-03;
  }
  if (r<82) {
    if (z<17)      return 0.034e-03;
    else if (z<70) return 0.07e-03;
    else           return 0.07e-03;
  }
  if (r<90) {
    if (z<17)      return 0.033e-03;
    else if (z<70) return 0.07e-03;
    else           return 0.07e-03;
  }
  if (r<100) {
    if (z<17)      return 0.033e-03;
    else if (z<70) return 0.07e-03;
    else           return 0.07e-03;
  }
  if (r<102) {
    if (z<17)      return 0.033e-03;
    else if (z<70) return 0.07e-03;
    else           return 0.07e-03;
  }
  return 0.;
}

#endif
