#ifndef _simulation_
#define _simulation_

#include "Track.h"
#include "Matrix.h"
#include "Propagation.h"

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, std::vector<Hit>& hits, int& charge, float pt);

float convertXYtoPhi(float X, float Y);
float convertXYtoPhiErr(float X, float Y, float XYerr);
float convertPhitoX(float phi, float X, float Y);
float convertPhitoY(float phi, float X, float Y);
SMatrixSym33 calcSimCov(float hitX, float hitY, float hitZ, float hitPosErrXY, float hitPosErrZ);

#endif
