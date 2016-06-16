#ifndef _PROPAGATION_KERNELS_H_
#define _PROPAGATION_KERNELS_H_

#include "GPlex.h"

void propagation_wrapper(cudaStream_t& stream,
    GPlexHV& msPar,
    GPlexLV& inPar, GPlexQI& inChg,
    GPlexLV& outPar, GPlexLL& errorProp,
    GPlexLS& outErr, 
    const int N);

void propagationForBuilding_wrapper(cudaStream_t& stream,
    const GPlexLS& inErr, const GPlexLV& inPar,
    const GPlexQI& inChg, const float radius,
    GPlexLS& outErr, GPlexLV& outPar, 
    const int N);

#endif  // _PROPAGATION_KERNELS_H_
