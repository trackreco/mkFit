#ifndef _conformalutils_
#define _conformalutils_

#include "Hit.h"
#include "Track.h"
#include "Matrix.h"

namespace mkfit {

void conformalFit(const Hit& hit0, const Hit& hit1, const Hit& hit2, TrackState& fitStateHit0, bool fiterrs = 1);

} // end namespace mkfit
#endif
