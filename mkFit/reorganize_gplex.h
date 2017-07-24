#ifndef REORGANIZE_GPLEX_H
#define REORGANIZE_GPLEX_H

#include "GPlex.h"
#include "Hit.h"
#include "HitStructuresCU.h"
#include "Track.h"

__device__ float *get_posArray(Hit &hit);
__device__ float *get_errArray(Hit &hit);
__device__ float *get_posArray(Track &track);
__device__ float *get_errArray(Track &track);

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

__device__ void OutputParErrCU_fn(Track *tracks, 
                                  const GPlexLS &Err, const GPlexLV &Par,
                                  const int beg, const int end, 
                                  const int itrack, const int N);

__device__ void OutputTracksCU_fn(Track *tracks, 
                                  const GPlexLS &Err_iP, const GPlexLV &Par_iP,
                                  const GPlexQI &Chg, const GPlexQF &Chi2,
                                  const GPlexQI &Label, const GPlexQI *HitsIdx,
                                  const int beg, const int end, 
                                  const int itrack, const int N,
                                  const bool update_hit_idx=true);

void InputTracksCU_wrapper(const cudaStream_t &stream, 
                           const EtaBinOfCandidatesCU &etaBin,
                           GPlexLS &Err_iP, GPlexLV &Par_iP,
                           GPlexQI &Chg, GPlexQF &Chi2, GPlexQI &Label,
                           GPlexQI *HitsIdx,
                           const int beg, const int end, const bool inputProp, int N);

void InputTracksAndHitsCU_wrapper(const cudaStream_t &stream, 
                                  Track *tracks, EventOfHitsCU &event_of_hits,
                                  GPlexLS &Err_iP, GPlexLV &Par_iP,
                                  GPlexHS *msErr_arr, GPlexHV *msPar_arr,
                                  GPlexQI &Chg, GPlexQF &Chi2, GPlexQI &Label,
                                  GPlexQI *HitsIdx,
                                  const int beg, const int end,
                                  const bool inputProp, int N);

void InputTracksAndHitIdxComb_wrapper(const cudaStream_t &stream, 
                                      const EtaBinOfCombCandidatesCU &etaBin,
                                      GPlexLS &Err_iP, GPlexLV &Par_iP,
                                      GPlexQI &Chg, GPlexQF &Chi2,
                                      GPlexQI &Label, GPlexQI *HitsIdx,
                                      GPlexQI &SeedIdx, GPlexQI &CandIdx,
                                      GPlexQB &Valid, const int Nhits,
                                      const int beg, const int end, 
                                      const bool inputProp, int N);

void OutputTracksCU_wrapper(const cudaStream_t &stream,
                            EtaBinOfCandidatesCU &etaBin,
                            GPlexLS &Err_iP, GPlexLV &Par_iP,
                            GPlexQI &Chg, GPlexQF &Chi2, GPlexQI &Label,
                            GPlexQI *HitsIdx,
                            const int beg, const int end, bool outputProp, int N);


void OutputFittedTracksCU_wrapper(const cudaStream_t &stream,
                                  Track *tracks_cu, 
                                  GPlexLS &Err_iP, GPlexLV &Par_iP,
                                  GPlexQI &Chg, GPlexQF &Chi2, GPlexQI &Label,
                                  const int beg, const int end, int N);

__device__
void GetHitErr(GPlexHS& msErr, const char* array, const int beg, const int end);
__device__
void GetHitPar(GPlexHV& msPar, const char* array, const int beg, const int end);


__device__ void InputTracksAndHitIdxComb_fn(Track *tracks, int *m_tracks_per_seed,
                                            GPlexLS &Err_iP, GPlexLV &Par_iP,
                                            GPlexQI &Chg, GPlexQF &Chi2,
                                            GPlexQI &Label, GPlexQI *HitsIdx,
                                            GPlexQI &SeedIdx, GPlexQI &CandIdx,
                                            GPlexQB &Valid,
                                            const int Nhits,
                                            const int beg, const int end, 
                                            const int itrack, const int N);

__device__ void InputTracksAndHitIdxComb_fn_seed(Track *tracks, int *m_tracks_per_seed,
                                            GPlexLS &Err_iP, GPlexLV &Par_iP,
                                            GPlexQI &Chg, GPlexQF &Chi2,
                                            GPlexQI &Label, GPlexQI *HitsIdx,
                                            GPlexQI &SeedIdx, GPlexQI &CandIdx,
                                            GPlexQB &Valid,
                                            const int Nhits,
                                            const int iseed_ev, 
                                            const int icand_ev,
                                            const int iseed_plex,
                                            const int N);

__device__ void OutputParErrCU_fn_seed(Track *tracks, 
                                  const GPlexLS &Err, const GPlexLV &Par,
                                  const int iseed_ev,
                                  const int icand_ev,
                                  int N);

#endif  // REORGANIZE_GPLEX_H
