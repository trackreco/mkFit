#ifndef _propagation_mplex_
#define _propagation_mplex_

#include "Track.h"
#include "Matrix.h"

void propagateLineToRMPlex(const MPlexLS &psErr,  const MPlexLV& psPar,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                                 MPlexLS &outErr,       MPlexLV& outPar,
                           const int      N_proc);

void propagateHelixToRMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const MPlexHV& msPar,
                                  MPlexLS &outErr,       MPlexLV& outPar,
                            const int      N_proc, const bool useParamBfield = false);

void propagateHelixToRMPlex(const MPlexLS& inErr,  const MPlexLV& inPar,
                            const MPlexQI& inChg,  const float    r,
                                  MPlexLS& outErr,       MPlexLV& outPar,
                            const int      N_proc);

void helixAtRFromIterativeCCSFullJac(const MPlexLV& inPar, const MPlexQI& inChg,
                                           MPlexLV& outPar, const MPlexQF &msRad,
                                           MPlexLL& errorProp,
                                     const int      N_proc);

void helixAtRFromIterativeCCS(const MPlexLV& inPar, const MPlexQI& inChg,
                                    MPlexLV& outPar, const MPlexQF &msRad,
                                    MPlexLL& errorProp,
                              const int      N_proc, const bool useParamBfield = false);

void helixAtRFromIterative(const MPlexLV& inPar, const MPlexQI& inChg, 
                                 MPlexLV& outPar, const MPlexQF &msRad, 
                                 MPlexLL& errorProp,
                           const int      N_proc, const bool useParamBfield = false);

void propagateHelixToZMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const MPlexHV& msPar,
                                  MPlexLS &outErr,       MPlexLV& outPar,
                            const int      N_proc, const bool useParamBfield = false);

void propagateHelixToZMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const float    z,
                                  MPlexLS &outErr,       MPlexLV& outPar,
                            const int      N_proc);

void helixAtZ(const MPlexLV& inPar,  const MPlexQI& inChg,
                    MPlexLV& outPar, const MPlexQF &msZ,
                    MPlexLL& errorProp,
              const int      N_proc, const bool useParamBfield = false);

void applyMaterialEffects(const MPlexQF &hitsRl, const MPlexQF& hitsXi, 
                                MPlexLS &outErr, MPlexLV& outPar,
                          const int      N_proc);

#endif
