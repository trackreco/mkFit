#ifndef _propagation_mplex_
#define _propagation_mplex_

#include "Track.h"
#include "Matrix.h"

void propagateLineToRMPlex(const MPlexLS &psErr,  const MPlexLV& psPar,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                                 MPlexLS &outErr,       MPlexLV& outPar);

void propagateHelixToRMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const MPlexHV& msPar,
			          MPlexLS &outErr,       MPlexLV& outPar);

void propagateHelixToRMPlex(const MPlexLS& inErr,  const MPlexLV& inPar,
                            const MPlexQI& inChg,  const float    r,
			          MPlexLS& outErr,       MPlexLV& outPar,
                            const int      N_proc);

//inline?
inline void computeJacobianSimple(int n, MPlexLL& errorProp, 
				  float s, float k, float p, float pxin, float pyin, float pzin, 
				  float TP, float cosTP, float sinTP);

void helixAtRFromIterative(const MPlexLV& inPar, const MPlexQI& inChg, 
			         MPlexLV& outPar, const MPlexQF &msRad, 
			         MPlexLL& errorProp, bool useSimpleJac);

void helixAtRFromIntersection(const MPlexLV& inPar, const MPlexQI& inChg, 
                                    MPlexLV& outPar, const MPlexQF &msRad, 
   			            MPlexLL& errorProp);

void applyMaterialEffects(const MPlexQF &hitsRl, const MPlexQF& hitsXi, 
			  MPlexLS &outErr, MPlexLV& outPar);

#endif
