#ifndef _PROPAGATION_KERNELS_H_
#define _PROPAGATION_KERNELS_H_

#include "GPlex.h"

void propagation_wrapper(cudaStream_t& stream,
    GPlex<float>& msPar,
    GPlex<float>& inPar, GPlex<int>& inChg,
    GPlex<float>& outPar, GPlex<float>& errorProp,
    GPlex<float>& outErr, 
    const int N);

#endif  // _PROPAGATION_KERNELS_H_
