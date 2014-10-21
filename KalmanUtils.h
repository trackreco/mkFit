#ifndef _kalmanutils_
#define _kalmanutils_

#include "Track.h"
#include "Matrix.h"

float computeChi2(TrackState& propagatedState, MeasurementState& measurementState, 
                  SMatrix36& projMatrix,SMatrix63& projMatrixT);

void zeroBlocksOutOf33(SMatrixSym66& matrix);
void copy33Into66(SMatrixSym33& in,SMatrixSym66& out);
void copy66Into33(SMatrixSym66& in,SMatrixSym33& out);
//see e.g. http://inspirehep.net/record/259509?ln=en
void updateParameters66(TrackState& propagatedState, MeasurementState& measurementState, TrackState& result);

//see e.g. http://inspirehep.net/record/259509?ln=en
TrackState updateParameters(TrackState& propagatedState, MeasurementState& measurementState, 
                            SMatrix36& projMatrix,SMatrix63& projMatrixT);

void updateParametersMPlex(const MPlexSS &psErr,  const MPlexMV& psPar,
                           const MPlexSS &msErr,  const MPlexMV& msPar,
                                 MPlexSS &outErr,       MPlexMV& outPar);

#endif
