#ifndef _buildtest_
#define _buildtest_

#include <map>
#include <string>
#include <vector>

#include "Track.h"
#include "Geometry.h"
#include "Event.h"

void buildTracksParallel(Event& ev);
void buildTracksSerial(Event& ev);
#endif
