#ifndef _fittest_mplex_
#define _fittest_mplex_

#include "Event.h"
#include "Track.h"

#ifdef USE_CUDA
#include "FitterCU.h"
#endif

double runFittingTestPlex(Event& ev, std::vector<Track>& rectracks);

#ifdef USE_CUDA
void runAllEventsFittingTestPlexGPU(std::vector<Event>& events);
double runFittingTestPlexGPU(FitterCU<float> &cuFitter, Event& ev, std::vector<Track>& rectracks);
#endif

#endif
