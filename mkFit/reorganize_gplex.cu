#include "reorganize_gplex.h"
#include <stdio.h>

#include "FitterCU.h"
#include "accessors_cu.h"
#include "Track.h"
#include "gpu_utils.h"

__device__ float *get_posArray(Hit &hit) {
    return hit.posArrayCU();
}
__device__ float *get_errArray(Hit &hit) {
    return hit.errArrayCU();
}

__device__ float *get_posArray(Track &track) {
    return track.posArrayCU();
}
__device__ float *get_errArray(Track &track) {
    return track.errArrayCU();
}

template <typename GPlexObj>
__device__ void SlurpIn_fn(GPlexObj to, // float *fArray, int stride, int kSize, 
                           const char *arr, const int *vi, const int N) {
  int j = threadIdx.x + blockDim.x * blockIdx.x;
  if (j<N) {
    const int *XHitPos = vi;
    const int off = XHitPos[j] * sizeof(Hit);
    for (int i = 0; i < to.kSize; ++i) { // plex_size
      /*fArray[i*stride+ j] = * (const T*) (arr + i*sizeof(T) + off);*/
      to(j, i, 0) = * (decltype(to.ptr)) (arr + i*sizeof(decltype(*to.ptr)) + off);
    }
  }
}


template <typename GPlexObj>
__device__ void SlurpInIdx_fn(GPlexObj to,
                             const char *arr, const int idx, const int N) {
  int j = threadIdx.x + blockDim.x * blockIdx.x;
  if (j<N) {
    for (int i = 0; i < to.kSize; ++i) { // plex_size
      to(j, i, 0) = * (decltype(to.ptr)) (arr + i*sizeof(decltype(*to.ptr)) + idx);
    }
  }
}


template <typename GPlexObj>
__device__ void SlurpOutIdx_fn(GPlexObj from, // float *fArray, int stride, int kSize, 
                               const char *arr, const int idx, const int N) {
  int j = threadIdx.x + blockDim.x * blockIdx.x;
  if (j<N) {
    for (int i = 0; i < from.kSize; ++i) { // plex_size
      * (decltype(from.ptr)) (arr + i*sizeof(decltype(*from.ptr)) + idx) = from(j, i , 0);
    }
  }
}


__device__ void HitToMs_fn(GPlexHS &msErr, GPlexHV &msPar,
                           Hit *hits, const GPlexQI &XHitSize, 
                           const GPlexHitIdx &XHitArr, 
                           GPlexQI &HitsIdx, const int hit_cnt, 
                           const int itrack, const int N) {
  /*int j = threadIdx.x + blockDim.x*blockIdx.x;*/
  //int itrack = threadIdx.x + blockDim.x * blockIdx.x;
  if (itrack < N) {

    const char *varr      = (char*) hits;
    const int   off_error = (char*) hits[0].errArrayCU() - varr;
    const int   off_param = (char*) hits[0].posArrayCU() - varr;

    if (hit_cnt < XHitSize[itrack]) {
      HitsIdx[itrack] = XHitArr(itrack, hit_cnt, 0) * sizeof(Hit);
    }
    SlurpInIdx_fn(msErr, varr + off_error, HitsIdx[itrack], N);
    SlurpInIdx_fn(msPar, varr + off_param, HitsIdx[itrack], N);
  }
}

__global__ void HitToMs_kernel(GPlexHS msErr, GPlexHV msPar, Hit *hits,
                               const GPlexQI XHitSize, const GPlexHitIdx XHitArr,
                               GPlexQI HitsIdx, const int hit_cnt, const int N) {
  int itrack = threadIdx.x + blockDim.x * blockIdx.x;
  HitToMs_fn(msErr, msPar, hits, XHitSize, XHitArr, HitsIdx, hit_cnt, itrack, N);
}

void HitToMs_wrapper(const cudaStream_t& stream,
                     GPlexHS &msErr, GPlexHV &msPar, LayerOfHitsCU &layer, 
                     const GPlexQI &XHitSize, const GPlexHitIdx &XHitArr,
                     GPlexQI &HitsIdx, int hit_cnt, const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  HitToMs_kernel <<< grid, block, 0 , stream >>>
    (msErr, msPar, layer.m_hits, XHitSize, XHitArr, HitsIdx, hit_cnt, N);
  cudaDeviceSynchronize();
}


__device__ void InputTracksCU_fn (Track *tracks, 
                                  GPlexLS &Err_iP, GPlexLV &Par_iP,
                                  GPlexQI &Chg, GPlexQF &Chi2,
                                  GPlexQI &Label, GPlexQI *HitsIdx,
                                  const int beg, const int end, 
                                  const int itrack, const int N) {
  //int itrack = threadIdx.x + blockDim.x*blockIdx.x;

  if (itrack < (end-beg) && itrack < N) {
    Track &trk = tracks[beg];
    const char *varr       = (char*) &trk;
    int   off_error = (char*) trk.errArrayCU() - varr;
    int   off_param = (char*) trk.posArrayCU() - varr;

    int i= itrack + beg;
    const Track &trk_i = tracks[i];
    int idx = (char*) &trk_i - varr;

    Label(itrack, 0, 0) = tracks[i].label();
    Chg(itrack, 0, 0) = tracks[i].charge();
    Chi2(itrack, 0, 0) = tracks[i].chi2();
    SlurpInIdx_fn(Err_iP, varr + off_error, idx, N);
    SlurpInIdx_fn(Par_iP, varr + off_param, idx, N);

    for (int hi = 0; hi < 3; ++hi)
      HitsIdx[hi](itrack, 0, 0) = tracks[i].getHitIdx(hi);//dummy value for now
  }
}


__global__ void InputTracksCU_kernel(Track *tracks, 
                                     GPlexLS Err_iP, GPlexLV Par_iP,
                                     GPlexQI Chg, GPlexQF Chi2, GPlexQI Label,
                                     GPlexQI *HitsIdx,
                                     int beg, int end, int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  InputTracksCU_fn(tracks, Err_iP, Par_iP, Chg, Chi2, Label, HitsIdx, beg, end, itrack, N);
}


void InputTracksCU_wrapper(const cudaStream_t &stream, 
                           const EtaBinOfCandidatesCU &etaBin,
                           GPlexLS &Err_iP, GPlexLV &Par_iP,
                           GPlexQI &Chg, GPlexQF &Chi2, GPlexQI &Label,
                           GPlexQI *HitsIdx,
                           const int beg, const int end, const bool inputProp, int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  InputTracksCU_kernel <<< grid, block, 0, stream >>>
    (etaBin.m_candidates, Err_iP, Par_iP, Chg, Chi2, Label, HitsIdx,
     beg, end, N);
}


__device__ void InputTracksAndHitsCU_fn (Track *tracks, LayerOfHitsCU *layerHits,
                                         GPlexLS &Err_iP, GPlexLV &Par_iP,
                                         GPlexHS *msErr_arr, GPlexHV *msPar_arr,
                                         GPlexQI &Chg, GPlexQF &Chi2,
                                         GPlexQI &Label, GPlexQI *HitsIdx,
                                         const int beg, const int end, 
                                         const int itrack, const int N) {
  //int itrack = threadIdx.x + blockDim.x*blockIdx.x;

  if (itrack < (end-beg) && itrack < N) {
    Track &trk = tracks[beg];
    const char *varr       = (char*) &trk;
    int   off_error = (char*) trk.errArrayCU() - varr;
    int   off_param = (char*) trk.posArrayCU() - varr;

    int i= itrack + beg;
    const Track &trk_i = tracks[i];
    int idx = (char*) &trk_i - varr;

    Label(itrack, 0, 0) = tracks[i].label();
    Chg(itrack, 0, 0) = tracks[i].charge();
    Chi2(itrack, 0, 0) = tracks[i].chi2();
    SlurpInIdx_fn(Err_iP, varr + off_error, idx, N);
    SlurpInIdx_fn(Par_iP, varr + off_param, idx, N);

    // Note Config::nLayers -- not suitable for building
    for (int hi = 0; hi < Config::nLayers; ++hi) {
      int hidx = tracks[i].getHitIdx(hi);
      Hit &hit = layerHits[hi].m_hits[hidx];

      HitsIdx[hi](itrack, 0, 0) = idx;
      if (hidx < 0) continue;

      SlurpInIdx_fn(msErr_arr[hi], (char *)hit.errArrayCU(), 0, N);
      SlurpInIdx_fn(msPar_arr[hi], (char *)hit.posArrayCU(), 0, N);
    }
  }
}


__global__ void InputTracksAndHitsCU_kernel(Track *tracks, LayerOfHitsCU *layers,
                                            GPlexLS Err_iP, GPlexLV Par_iP,
                                            GPlexHS *msErr_arr, GPlexHV *msPar_arr,
                                            GPlexQI Chg, GPlexQF Chi2, GPlexQI Label,
                                            GPlexQI *HitsIdx,
                                            int beg, int end, int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  InputTracksAndHitsCU_fn(tracks, layers, Err_iP, Par_iP, msErr_arr, msPar_arr,
                          Chg, Chi2, Label, HitsIdx, beg, end, itrack, N);
}


void InputTracksAndHitsCU_wrapper(const cudaStream_t &stream, 
                                  Track *tracks, EventOfHitsCU &event_of_hits,
                                  GPlexLS &Err_iP, GPlexLV &Par_iP,
                                  GPlexHS *msErr_arr, GPlexHV *msPar_arr,
                                  GPlexQI &Chg, GPlexQF &Chi2, GPlexQI &Label,
                                  GPlexQI *HitsIdx,
                                  const int beg, const int end, 
                                  const bool inputProp, int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  InputTracksAndHitsCU_kernel <<< grid, block, 0, stream >>>
    (tracks, event_of_hits.m_layers_of_hits, 
     Err_iP, Par_iP, 
     msErr_arr, msPar_arr, 
     Chg, Chi2, Label, HitsIdx,
     beg, end, N);
}


__device__ void OutputTracksCU_fn(Track *tracks, 
                                  const GPlexLS &Err_iP, const GPlexLV &Par_iP,
                                  const GPlexQI &Chg, const GPlexQF &Chi2,
                                  const GPlexQI &Label, const GPlexQI *HitsIdx,
                                  const int beg, const int end, 
                                  const int itrack, const int N,
                                  const bool update_hit_idx) {
  //int itrack = threadIdx.x + blockDim.x*blockIdx.x;

  if (itrack < (end-beg) && itrack < N) {
    Track &trk = tracks[beg];
    const char *varr       = (char*) &trk;
    int   off_error = (char*) trk.errArrayCU() - varr;
    int   off_param = (char*) trk.posArrayCU() - varr;

    int i= itrack + beg;
    const Track &trk_i = tracks[i];
    int idx = (char*) &trk_i - varr;

    SlurpOutIdx_fn(Err_iP, varr + off_error, idx, N);
    SlurpOutIdx_fn(Par_iP, varr + off_param, idx, N);
    tracks[i].setCharge(Chg(itrack, 0, 0));
    tracks[i].setChi2(Chi2(itrack, 0, 0));
    tracks[i].setLabel(Label(itrack, 0, 0));

    if (update_hit_idx) {
      tracks[i].resetHits();
      /*int nGoodItIdx = 0;*/
      for (int hi = 0; hi < Config::nLayers; ++hi) {
        tracks[i].addHitIdx(HitsIdx[hi](itrack, 0, 0),0.);
        // FIXME: We probably want to use registers instead of going for class members:
        /*int hit_idx = HitsIdx[hi](itrack, 0, 0);*/
        /*tracks[i].setHitIdx(hi, hit_idx);*/
        /*if (hit_idx >= 0) {*/
        /*nGoodItIdx++; */
        /*}*/
      }
      /*tracks[i].setNGoodHitIdx(nGoodItIdx);*/
      /*tracks[i].setChi2(0.);*/
    }
  }
}

__global__ void OutputTracksCU_kernel(Track *tracks, 
                                     GPlexLS Err_iP, GPlexLV Par_iP,
                                     GPlexQI Chg, GPlexQF Chi2, GPlexQI Label,
                                     GPlexQI *HitsIdx,
                                     int beg, int end, int N,
                                     const bool update_hit_idx=true) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  OutputTracksCU_fn(tracks, Err_iP, Par_iP, Chg, Chi2, Label, HitsIdx,
                    beg, end, itrack, N, update_hit_idx);
}


void OutputTracksCU_wrapper(const cudaStream_t &stream,
                            EtaBinOfCandidatesCU &etaBin,
                            GPlexLS &Err_iP, GPlexLV &Par_iP,
                            GPlexQI &Chg, GPlexQF &Chi2, GPlexQI &Label,
                            GPlexQI *HitsIdx,
                            const int beg, const int end, const bool outputProp, int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  OutputTracksCU_kernel <<< grid, block, 0, stream >>>
    (etaBin.m_candidates, Err_iP, Par_iP, Chg, Chi2, Label, HitsIdx, beg, end, N);
}


void OutputFittedTracksCU_wrapper(const cudaStream_t &stream,
                                  Track *tracks_cu, 
                                  GPlexLS &Err_iP, GPlexLV &Par_iP,
                                  GPlexQI &Chg, GPlexQF &Chi2, GPlexQI &Label,
                                  const int beg, const int end, int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  OutputTracksCU_kernel <<< grid, block, 0, stream >>>
    (tracks_cu, Err_iP, Par_iP, Chg, Chi2, Label, nullptr, beg, end, N, false);
}
