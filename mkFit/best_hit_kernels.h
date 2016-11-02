#ifndef BEST_HIT_KERNELS_H
#define BEST_HIT_KERNELS_H 


#include "GPlex.h"
#include "HitStructuresCU.h"
#include "GeometryCU.h"


void getNewBestHitChi2_wrapper(const cudaStream_t &stream,
    const GPlexQI &XHitSize, const GPlexHitIdx &XHitArr,
    const GPlexQF &outChi2, float *minChi2, int *bestHit,
    const int hit_cnt, const int N);


void fill_array_cu(float *array, const int size, const float value);


void updateTracksWithBestHit_wrapper(const cudaStream_t &stream,
    LayerOfHitsCU &, const float *minChi2, const int *best_hit, 
    GPlexHS &msErr, GPlexHV &msPar, const GPlexLV &propPar,
    GPlexQF &Chi2, GPlexQI& HitsIdx, const int N);

int getMaxNumHits_wrapper(const GPlexQI d_XHitSize, const int N);


void bestHit_wrapper(const cudaStream_t &stream,
    LayerOfHitsCU &layer, const GPlexQI &XHitSize, const GPlexHitIdx &XHitArr,
    const GPlexLS &propErr, GPlexHS &msErr, GPlexHV &msPar,
    const GPlexLV &propPar, GPlexQF &outChi2,
    GPlexQF &Chi2, GPlexQI& HitsIdx,
    const int maxSize, const int N);


void findBestHit_wrapper(cudaStream_t &stream,
                         LayerOfHitsCU *layers, 
                         EventOfCandidatesCU &event_of_cands_cu,
                         GPlexQI &XHitSize, GPlexHitIdx &XHitArr,
                         GPlexLS &Err_iP, GPlexLV &Par_iP,
                         GPlexHS *msErr, GPlexHV *msPar,
                         GPlexLS &Err_iC, GPlexLV &Par_iC,
                         GPlexQF &outChi2,
                         GPlexQF &Chi2, GPlexQI *HitsIdx,
                         GPlexQI &inChg, GPlexQI &Label,
                         GeometryCU &geom, 
                         int *maxSize, int N);


#endif  /* ifndef BEST_HIT_KERNELS_H */
