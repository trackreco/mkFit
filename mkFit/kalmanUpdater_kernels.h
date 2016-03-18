#ifndef _KALMAN_UPDATER_KERNELS_H_
#define _KALMAN_UPDATER_KERNELS_H_

#include "GPlex.h"

void kalmanUpdate_wrapper(cudaStream_t& stream,
    GPlex<float>& d_propErr, GPlex<float>& d_msErr,
    GPlex<float>& d_par_iP, GPlex<float>& d_msPar,
    GPlex<float>& d_par_iC, GPlex<float>& d_outErr,
    const int N);

void reorganizeMs_wrapper(cudaStream_t& stream, GPlex<float>& msPar, float *full_posArray,
    GPlex<float>& msErr, float *full_errArray, int *full_hitIdx, int hi, int maxHits,
    int N, int hs, int hv, int Nhits);

#endif  // _KALMAN_UPDATER_KERNELS_H_
