#ifndef _KALMAN_UPDATER_KERNELS_H_
#define _KALMAN_UPDATER_KERNELS_H_

#include "GPlex.h"

#if 0
void addIntoUpperLeft3x3_wrapper(dim3 grid, dim3 block, cudaStream_t& stream,
    GPlex<float>& d_propErr,
    GPlex<float>& msErr, GPlex<float>& d_resErr, const int N);

void invertCramerSym_wrapper(dim3 grid, dim3 block, cudaStream_t& stream, GPlex<float>& d_resErr, const int N);

void upParam_MultKalmanGain_wrapper(dim3 grid, dim3 block, cudaStream_t& stream,
    GPlex<float>& d_Err_iP, GPlex<float>& d_resErr, GPlex<float>& d_kalmanGain, int num_matrices);

void multResidualsAdd_wrapper(dim3 grid, dim3 block, cudaStream_t& stream,
    GPlex<float>& d_kalmanGain, GPlex<float>& d_par_iP,
    GPlex<float>& d_msPar, GPlex<float>& d_par_iC, int N);

void kalmanGain_x_propErr_wrapper(dim3 grid, dim3 block, cudaStream_t& stream,
    GPlex<float>& d_kalmanGain, GPlex<float>& d_propErr, 
    GPlex<float>& d_outErr, const int num_matrices);
#endif

// TODO: temporary function. merge routines to reduce the number of global mem. accesses.
void kalmanUpdateMerged_wrapper(cudaStream_t& stream,
    GPlex<float>& d_propErr, GPlex<float>& d_msErr,
    GPlex<float>& d_kalmanGain,
    GPlex<float>& d_par_iP, GPlex<float>& d_msPar,
    GPlex<float>& d_par_iC, GPlex<float>& d_outErr,
    const int N);

void reorganizeMs_wrapper(cudaStream_t& stream, GPlex<float>& msPar, float *full_posArray,
    GPlex<float>& msErr, float *full_errArray, int *full_hitIdx, int hi, int maxHits,
    int N, int hs, int hv, int Nhits);
#endif  // _KALMAN_UPDATER_KERNELS_H_
