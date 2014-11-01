#ifndef _simulation_
#define _simulation_

#include "Track.h"
#include "Matrix.h"
#include "Propagation.h"
#include "Geometry.h"

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, HitVec& hits, int& charge, float pt, Geometry*, HitVec* initHits = nullptr);
#endif
