#ifndef _kalmanutils_
#define _kalmanutils_

#include "Track.h"
#include "Matrix.h"

//float computeChi2(const TrackState& propagatedState, const MeasurementState& measurementState);
inline float computeChi2(const TrackState& propagatedState, const MeasurementState& measurementState) {
  /*
  int invFail(0);
  const SMatrix33 resErr = measurementState.errors() + propagatedState.errors.Sub<SMatrix33>(0,0);
  return ROOT::Math::Similarity(measurementState.parameters()-propagatedState.parameters.Sub<SVector3>(0),
                                resErr.InverseFast(invFail));
  */
  float r = getHypot(measurementState.pos_[0],measurementState.pos_[1]);
  //rotate to the tangent plane to the cylinder of radius r at the hit position
  SMatrix33 rot;
  rot[0][0] = -measurementState.pos_[1]/r;
  rot[0][1] = 0;
  rot[0][2] = measurementState.pos_[0]/r;
  rot[1][0] = rot[0][2];
  rot[1][1] = 0;
  rot[1][2] = -rot[0][0];
  rot[2][0] = 0;
  rot[2][1] = 1;
  rot[2][2] = 0;
  const SMatrix33 rotT = ROOT::Math::Transpose(rot);
  const SVector3 res_glo = measurementState.parameters()-propagatedState.parameters.Sub<SVector3>(0);
  const SVector3 res_loc3 = rotT * res_glo;
  //the matrix to invert has to be 2x2
  const SVector2 res(res_loc3[0],res_loc3[1]);
  const SMatrixSym33 resErr_glo = measurementState.errors() + propagatedState.errors.Sub<SMatrixSym33>(0,0);
  const SMatrixSym22 resErr = ROOT::Math::SimilarityT(rot,resErr_glo).Sub<SMatrixSym22>(0,0);
  int invFail(0);
  SMatrixSym22 resErrInv = resErr.InverseFast(invFail);
  return ROOT::Math::Similarity(res,resErrInv);
}

//see e.g. http://inspirehep.net/record/259509?ln=en
void updateParameters66(TrackState& propagatedState, MeasurementState& measurementState, TrackState& result);
TrackState updateParameters(const TrackState& propagatedState, const MeasurementState& measurementState);

#endif
