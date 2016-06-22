#ifndef _kalmanutils_mplex_
#define _kalmanutils_mplex_

#include "Track.h"
#include "Matrix.h"

void updateParametersMPlex(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                                 MPlexLS &outErr,       MPlexLV& outPar,
                           const int      N_proc);

void computeChi2MPlex(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
		      const MPlexHS &msErr,  const MPlexHV& msPar,
                            MPlexQF& outChi2,
                      const int      N_proc);

void updateParametersEndcapMPlex(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
				 const MPlexHS &msErr,  const MPlexHV& msPar,
                                       MPlexLS &outErr,       MPlexLV& outPar,
				 const int      N_proc);

void computeChi2EndcapMPlex(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
			    const MPlexHS &msErr,  const MPlexHV& msPar,
                                  MPlexQF& outChi2,
			    const int      N_proc);

#endif
