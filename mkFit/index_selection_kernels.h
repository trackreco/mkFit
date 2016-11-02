#ifndef _INDEX_SELECTION_KERNELS_H_
#define _INDEX_SELECTION_KERNELS_H_

#include "HitStructuresCU.h"
#include "GPlex.h"

void selectHitIndices_wrapper(const cudaStream_t& stream,
    const LayerOfHitsCU& layer_of_hits, const GPlexLS& Err, const GPlexLV& Par, 
    GPlexQI& XHitSize, GPlexHitIdx& XHitArr, const int N);

__device__ void selectHitIndices_fn(const LayerOfHitsCU &layer_of_hits,
    const GPlexLS &Err, const GPlexLV &Par, GPlexQI &XHitSize,
    GPlexHitIdx &XHitArr, const int itrack, const int N);

#endif  // _INDEX_SELECTION_KERNELS_H_
