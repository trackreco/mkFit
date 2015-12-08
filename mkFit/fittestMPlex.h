#ifndef _fittest_mplex_
#define _fittest_mplex_

#include "Event.h"
#include "Track.h"

void   make_validation_tree(const char         *fname,
                            std::vector<Track> &simtracks,
                            std::vector<Track> &rectracks);

double runFittingTestPlex(Event& ev, std::vector<Track>& rectracks);

double runFittingTestPlexGPU(Event& ev, std::vector<Track>& rectracks);
#endif
