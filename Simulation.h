#ifndef _simulation_
#define _simulation_

#include "Track.h"
#include "Matrix.h"
#include "Propagation.h"
#include "Geometry.h"

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, std::vector<Hit>& hits, int& charge, float pt, Geometry*, std::vector<Hit>& initHits);
void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, std::vector<Hit>& hits, int& charge, float pt, Geometry*);

#endif
