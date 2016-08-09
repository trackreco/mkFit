#ifndef REORGANIZE_GPLEX_H
#define REORGANIZE_GPLEX_H

#include "GPlex.h"
#include "Hit.h"
#include "HitStructuresCU.h"

__device__ void HitToMs_fn(GPlexHS &msErr, GPlexHV &msPar,
                           Hit *hits, const GPlexQI &XHitSize,
                           const GPlexHitIdx &XHitArr, 
                           GPlexQI &HitsIdx, const int hit_cnt,
                           const int itrack, const int N);

__global__ void HitToMs_kernel(GPlexHS msErr, GPlexHV msPar, Hit *hits, 
                               const GPlexQI XHitSize, const GPlexHitIdx XHitArr, 
                               GPlexQI HitsIdx, const int hit_cnt, const int N);

void HitToMs_wrapper(const cudaStream_t& stream,
                     GPlexHS &msErr, GPlexHV &msPar, LayerOfHitsCU &layer, 
                     const GPlexQI &XHitSize, const GPlexHitIdx &XHitArr, 
                     GPlexQI &HitsIdx, const int hit_cnt, const int N);

__device__ void InputTracksCU_fn(Track *tracks, 
                                 GPlexLS &Err_iP, GPlexLV &Par_iP,
                                 GPlexQI &Chg, GPlexQF &Chi2,
                                 GPlexQI &Label, GPlexQI *HitsIdx,
                                 const int beg, const int end,
                                 const int itrack, const int N);

__device__ void OutputTracksCU_fn(Track *tracks, 
                                  const GPlexLS &Err_iP, const GPlexLV &Par_iP,
                                  const GPlexQI &Chg, const GPlexQF &Chi2,
                                  const GPlexQI &Label, const GPlexQI *HitsIdx,
                                  const int beg, const int end, 
                                  const int itrack, const int N);

void InputTracksCU_wrapper(const cudaStream_t &stream, 
                           const EtaBinOfCandidatesCU &etaBin,
                           GPlexLS &Err_iP, GPlexLV &Par_iP,
                           GPlexQI &Chg, GPlexQF &Chi2, GPlexQI &Label,
                           GPlexQI *HitsIdx,
                           const int beg, const int end, const bool inputProp, int N);

void OutputTracksCU_wrapper(const cudaStream_t &stream,
                            EtaBinOfCandidatesCU &etaBin,
                            GPlexLS &Err_iP, GPlexLV &Par_iP,
                            GPlexQI &Chg, GPlexQF &Chi2, GPlexQI &Label,
                            GPlexQI *HitsIdx,
                            const int beg, const int end, const bool outputProp, int N);

#endif  // REORGANIZE_GPLEX_H
