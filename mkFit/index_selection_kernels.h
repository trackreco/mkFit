#ifndef _INDEX_SELECTION_KERNELS_H_
#define _INDEX_SELECTION_KERNELS_H_

#include "HitStructuresCU.h"
#include "GPlex.h"

void selectHitIndices_wrapper(cudaStream_t& stream,
    LayerOfHitsCU& layer_of_hits, GPlexLS& Err, GPlexLV& Par, 
    GPlexQI& XHitSize, GPlexHitIdx& XHitArr, int N);

#endif  // _INDEX_SELECTION_KERNELS_H_
