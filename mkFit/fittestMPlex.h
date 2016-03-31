#ifndef _fittest_mplex_
#define _fittest_mplex_

#include "Event.h"
#include "Track.h"

#ifdef USE_CUDA
#include "FitterCU.h"
#endif

void   make_validation_tree(const char         *fname,
                            std::vector<Track> &simtracks,
                            std::vector<Track> &rectracks);

double runFittingTestPlex(Event& ev, std::vector<Track>& rectracks);

#ifdef USE_CUDA
void runAllEventsFittingTestPlexGPU(std::vector<Event>& events);
double runFittingTestPlexGPU(FitterCU<float> &cuFitter, Event& ev, std::vector<Track>& rectracks);
#endif

#endif
