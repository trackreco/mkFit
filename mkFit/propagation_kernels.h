#ifndef _PROPAGATION_KERNELS_H_
#define _PROPAGATION_KERNELS_H_

#include "GPlex.h"

void helixAtRFromIterative_wrapper(dim3 grid, dim3 block, cudaStream_t& stream, 
    GPlex<float>& inPar, GPlex<int>& inChg, GPlex<float>& outPar, 
    GPlex<float>& msRad, GPlex<float>& errorProp, int N);

void similarity_wrapper(dim3 grid, dim3 block, cudaStream_t& stream,
    GPlex<float>& errorProp, GPlex<float>& outErr, int N);

void computeMsRad_wrapper(dim3 grid, dim3 block, cudaStream_t& stream,
    GPlex<float>& msPar, GPlex<float>& msRad, int cN);

void propagationMerged_wrapper(cudaStream_t& stream,
    GPlex<float>& msPar,
    GPlex<float>& inPar, GPlex<int>& inChg,
    GPlex<float>& outPar, GPlex<float>& errorProp,
    GPlex<float>& outErr, 
    const int N);

#endif  // _PROPAGATION_KERNELS_H_
