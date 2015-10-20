#ifndef _kalmanutils_
#define _kalmanutils_

#include "Track.h"
#include "Matrix.h"

//float computeChi2(const TrackState& propagatedState, const MeasurementState& measurementState);
inline float computeChi2(const TrackState& propagatedState, const MeasurementState& measurementState) {
  int invFail(0);
  const SMatrix33 resErr = measurementState.errors() + propagatedState.errors.Sub<SMatrix33>(0,0);
  return ROOT::Math::Similarity(measurementState.parameters()-propagatedState.parameters.Sub<SVector3>(0),
                                resErr.InverseFast(invFail));
}

//see e.g. http://inspirehep.net/record/259509?ln=en
void updateParameters66(TrackState& propagatedState, MeasurementState& measurementState, TrackState& result);
TrackState updateParameters(const TrackState& propagatedState, const MeasurementState& measurementState);

#endif
