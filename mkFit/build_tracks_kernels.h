#ifndef _BUILD_TRACKS_KERNELS_
#define _BUILD_TRACKS_KERNELS_

#include "GPlex.h"
#include "HitStructuresCU.h"

void getHitFromLayer_wrappper( const cudaStream_t& stream,
    LayerOfHitsCU& layer_cu, GPlexQI& HitsIdx, 
    GPlexHV& msPar, GPlexHS& msErr, int N);
void UpdateMissingHits_wrapper(
    const cudaStream_t& stream, GPlexQI& HitsIdx, 
    GPlexLV& Par_iP, GPlexLS& Err_iP,
    GPlexLV& Par_iC, GPlexLS& Err_iC,
    int N);

__device__
void UpdateWithLastHit_fn(
    LayerOfHitsCU& layer_of_hits, GPlexQI& HitsIdx, 
    GPlexHV& msPar,  GPlexHS& msErr,
    GPlexLV& Par_iP, GPlexLS& Err_iP,
    GPlexLV& Par_iC, GPlexLS& Err_iC,
    int itrack_plex, int N);

void UpdateWithLastHit_wrapper(
    const cudaStream_t& stream,
    LayerOfHitsCU& layer_cu, GPlexQI& HitsIdx, 
    GPlexHV& msPar, GPlexHS& msErr,
    GPlexLV& Par_iP, GPlexLS& Err_iP,
    GPlexLV& Par_iC, GPlexLS& Err_iC,
    int N);

#endif  /* ifndef _BUILD_TRACKS_KERNELS_ */ 
