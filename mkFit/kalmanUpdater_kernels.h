#ifndef _KALMAN_UPDATER_KERNELS_H_
#define _KALMAN_UPDATER_KERNELS_H_

#include "GPlex.h"

void kalmanUpdate_wrapper(cudaStream_t& stream,
    GPlexLS& d_propErr, GPlexHS& d_msErr,
    GPlexLV& d_par_iP, GPlexHV& d_msPar,
    GPlexLV& d_par_iC, GPlexLS& d_outErr,
    const int N);

void reorganizeMs_wrapper(cudaStream_t& stream, GPlexQF& msPar,
    float *full_posArray, GPlexHS& msErr, 
    float *full_errArray, int *full_hitIdx, int hi, int maxHits,
    int N, int hs, int hv, int Nhits);

__device__ void addIntoUpperLeft3x3_fn(const GPlexLS __restrict__ &A,
                                       const GPlexHS __restrict__ &B,
                                       GPlexRegHS &c, const int N, int n);

__device__ void subtractFirst3_fn(const GPlexHV __restrict__ &A,
                                  const GPlexLV __restrict__ &B,
                                  GPlexRegHV &C, const int N, int n);

__device__ void invertCramerSym_fn(float *a);
__device__ void invertCramerSym2x2_fn(GPlexReg2S &a);

#endif  // _KALMAN_UPDATER_KERNELS_H_
