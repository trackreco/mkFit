#ifndef _kalmanutils_mplex_
#define _kalmanutils_mplex_

#include "Track.h"
#include "Matrix.h"

#ifndef __APPLE__

void updateParametersMPlex(const MPlexLS &psErr,  const MPlexLV& psPar,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                                 MPlexLS &outErr,       MPlexLV& outPar);
#endif

#endif
