#ifndef _fittest_mplex_
#define _fittest_mplex_

#include "Track.h"

void   generateTracks(std::vector<Track>& simtracks, int Ntracks);
void   make_validation_tree(const char         *fname,
                            std::vector<Track> &simtracks,
                            std::vector<Track> &rectracks);

double runFittingTest(std::vector<Track>& simtracks, std::vector<Track>& rectracks);

#ifndef __APPLE__
double runFittingTestPlex(std::vector<Track>& simtracks, std::vector<Track>& rectracks);
#endif

#endif
