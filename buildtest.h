#ifndef _buildtest_
#define _buildtest_

#include <map>
#include <string>
#include <vector>

#include "Track.h"
#include "Geometry.h"
#include "Event.h"
#include "BinInfoUtils.h"

void buildTracksBySeeds(const BinInfoMap &, Event& ev);
void buildTracksByLayers(const BinInfoMap &, Event& ev);
#endif
