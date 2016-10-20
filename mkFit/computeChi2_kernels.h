#ifndef _COMPUTE_CHI2_KERNELS_H_
#define _COMPUTE_CHI2_KERNELS_H_

#include "HitStructuresCU.h"
#include "GPlex.h"
#include "GeometryCU.h"


__device__ void computeChi2_fn(const GPlexLS &propErr, const GPlexHS &msErr,
    const GPlexHV &msPar, const GPlexLV &propPar, GPlexQF &outChi2, 
    const int itrack, const int N);


void computeChi2_wrapper(cudaStream_t &stream, 
    const GPlexLS &propErr, const GPlexHS &msErr, const GPlexHV &msPar, 
    const GPlexLV &propPar, GPlexQF &outChi2, const int N);

#if 1
void HitToMs_wrapper(cudaStream_t& stream,
    GPlexHS &msErr, GPlexHV &msPar, LayerOfHitsCU &layer, GPlexQI &XHitSize, 
    GPlexHitIdx &XHitArr, GPlexQI &HitsIdx, int hit_cnt, int N);
#endif

__device__ void RotateResidulsOnTangentPlane_fn(const GPlexRegQF& r00, 
    const GPlexRegQF& r01, const GPlexRegHV &a, GPlexReg2V &b);

__device__ void ProjectResErr_fn(const GPlexRegQF& a00, const GPlexRegQF& a01,
                                 const GPlexRegHS &b, GPlexRegHH &c);

__device__ void ProjectResErrTransp_fn(const GPlexRegQF& a00, const GPlexRegQF& a01, 
                                       const GPlexRegHH &b, GPlexReg2S &c);

#endif
