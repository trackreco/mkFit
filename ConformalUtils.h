#ifndef _conformalutils_
#define _conformalutils_

#include "Track.h"
#include "Matrix.h"

void conformalFit(const Hit& hit0, const Hit& hit1, const Hit& hit2, int charge, TrackState& fitStateHit0, bool fiterrs = 1);

#endif
