#ifndef _buildtest_
#define _buildtest_

#include <map>
#include <string>
#include <vector>

#include "Track.h"
#include "Geometry.h"
#include "Event.h"

void buildTracks(Event& ev, const int nlayers_per_seed,
  const unsigned int maxCand, const float chi2Cut, const float nSigma, const float minDPhi);
#endif
