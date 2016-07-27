#ifndef _kalmanutils_
#define _kalmanutils_

#include "Track.h"
#include "Matrix.h"

float computeChi2(const TrackState& propagatedState, const MeasurementState& measurementState);

//see e.g. http://inspirehep.net/record/259509?ln=en
void updateParameters66(TrackState& propagatedState, MeasurementState& measurementState, TrackState& result);
TrackState updateParameters(const TrackState& propagatedState, const MeasurementState& measurementState);

TrackState updateParametersEndcap(const TrackState& propagatedState, const MeasurementState& measurementState);
float computeChi2Endcap(const TrackState& propagatedState, const MeasurementState& measurementState);

#endif
