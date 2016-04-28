#ifndef _conformalutils_mplex_
#define _conformalutils_mplex_

#include "Hit.h"
#include "Track.h"
#include "Matrix.h"

  // write to iC --> next step will be a propagation no matter what
void conformalFitMPlex(bool fitting, const MPlexQI inChg, 
		       MPlexLS& outErr, MPlexLV& outPar, 
		       const MPlexHV& msPar0, const MPlexHV& msPar1, const MPlexHV& msPar2);

#endif
