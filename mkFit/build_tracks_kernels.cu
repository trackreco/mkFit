#include "build_tracks_kernels.h"

#include "GPlex.h"
#include "reorganize_gplex.h"
#include "kalmanUpdater_kernels.h"
#include "computeChi2_kernels.h"

constexpr int BLOCK_SIZE_X = 16;


__device__ void getHitFromLayer_fn(LayerOfHitsCU& layer_of_hits, 
    GPlexQI& HitsIdx, GPlexHV& msPar, GPlexHS& msErr, int itrack_plex, int N)
{
  if (itrack_plex < N)
  {
    int hit_idx = HitsIdx[itrack_plex];
    if (hit_idx >= 0)
    {
      Hit &hit = layer_of_hits.m_hits[hit_idx];

      GetHitErr(msErr, (char *)hit.errArrayCU(), 0, N);
      GetHitPar(msPar, (char *)hit.posArrayCU(), 0, N);
    }
  }
}

__global__ void getHitFromLayer_kernel(LayerOfHitsCU& layer_of_hits, 
    GPlexQI HitsIdx, GPlexHV msPar, GPlexHS msErr, int N)
{
  int itrack_plex = threadIdx.x + blockDim.x * blockIdx.x;
  getHitFromLayer_fn(layer_of_hits, HitsIdx, msPar, msErr, itrack_plex, N);
}

void getHitFromLayer_wrappper( const cudaStream_t& stream,
    LayerOfHitsCU& layer_cu, GPlexQI& HitsIdx, 
    GPlexHV& msPar, GPlexHS& msErr, int N)
{
  int gridx = (N-1)/BLOCK_SIZE_X + 1;
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  getHitFromLayer_kernel <<< grid, block, 0, stream >>>
    (layer_cu, HitsIdx, msPar, msErr, N);
}


__device__ void updateMissingHits_fn(GPlexQI& HitsIdx, 
    GPlexLV& Par_iP, GPlexLS& Err_iP,
    GPlexLV& Par_iC, GPlexLS& Err_iC, int i, int N)
{
  if (i < N)
  {
    if (HitsIdx[i] < 0)
    {
      for (auto j = 0; j < Err_iC.kSize; ++j) {
        Err_iC[i + j * Err_iC.stride] = Err_iP[i + j * Err_iP.stride];
      }
      for (auto j = 0; j < Par_iC.kSize; ++j) {
        Par_iC[i + j * Par_iC.stride] = Par_iP[i + j * Par_iP.stride];
      }
    }
  }
}


__global__ void updateMissingHits_kernel(GPlexQI HitsIdx, 
    GPlexLV Par_iP, GPlexLS Err_iP,
    GPlexLV Par_iC, GPlexLS Err_iC, int N)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  updateMissingHits_fn(HitsIdx, Par_iP, Err_iP, Par_iC, Err_iC, i, N);
}


void UpdateMissingHits_wrapper(
    const cudaStream_t& stream, GPlexQI& HitsIdx, 
    GPlexLV& Par_iP, GPlexLS& Err_iP,
    GPlexLV& Par_iC, GPlexLS& Err_iC,
    int N)
{
  int gridx = (N-1)/BLOCK_SIZE_X + 1;
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  updateMissingHits_kernel <<< grid, block, 0, stream >>>
    (HitsIdx, Par_iP, Err_iP, Par_iC, Err_iC, N);
}


__device__
void UpdateWithLastHit_fn(
    LayerOfHitsCU& layer_of_hits, GPlexQI& HitsIdx, 
    GPlexHV& msPar,  GPlexHS& msErr,
    GPlexLV& Par_iP, GPlexLS& Err_iP,
    GPlexLV& Par_iC, GPlexLS& Err_iC,
    int itrack_plex, int N)
{
  getHitFromLayer_fn(layer_of_hits, HitsIdx, msPar, msErr, itrack_plex, N);
  kalmanUpdate_fn(Err_iP, msErr, Par_iP, msPar, Par_iC, Err_iC, itrack_plex, N);
  updateMissingHits_fn(HitsIdx, Par_iP, Err_iP, Par_iC, Err_iC, itrack_plex, N);
}


__global__
void UpdateWithLastHit_kernel(
    LayerOfHitsCU& layer_of_hits, GPlexQI HitsIdx, 
    GPlexHV msPar, GPlexHS msErr,
    GPlexLV Par_iP, GPlexLS Err_iP,
    GPlexLV Par_iC, GPlexLS Err_iC,
    int N)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  UpdateWithLastHit_fn(layer_of_hits, HitsIdx, msPar, msErr, 
                       Par_iP, Err_iP, Par_iC, Err_iC, i, N);
}


void UpdateWithLastHit_wrapper(
    const cudaStream_t& stream,
    LayerOfHitsCU& layer_cu, GPlexQI& HitsIdx, 
    GPlexHV& msPar, GPlexHS& msErr,
    GPlexLV& Par_iP, GPlexLS& Err_iP,
    GPlexLV& Par_iC, GPlexLS& Err_iC,
    int N)
{
  int gridx = (N-1)/BLOCK_SIZE_X + 1;
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  UpdateWithLastHit_kernel <<< grid, block, 0, stream >>>
    (layer_cu, HitsIdx, msPar, msErr,
     Par_iP, Err_iP, Par_iC, Err_iC, N);
}


#if 0
__device__ void findCandidates_fn(
    Hit *hits, const GPlexQI &XHitSize, const GPlexHitIdx &XHitArr, 
    const GPlexLS &propErr, GPlexHS &msErr, GPlexHV &msPar,
    const GPlexLV &propPar, GPlexQF &outChi2,
    GPlexQF &Chi2, GPlexQI &HitsIdx,
    const int maxSize, const int itrack, const int N) {

  // FIXME: We probably don't want to find candidates for
  //        non-existing tracks. 
  //        However, Valid will need to be updated at some point.
  //    This is just a note to remember to do it.
  //if (!Valid[itrack]) return;

  if (itrack < N)
    HitsIdx[itrack] = 0;

  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
    HitToMs_fn(msErr, msPar, hits, XHitSize, XHitArr, HitsIdx, hit_cnt, itrack, N);
    computeChi2_fn(propErr, msErr, msPar, propPar, outChi2, itrack, N);
    // get the maxCandsPerSeed best cands for this TRACK
    //getNewBestHitChi2_fn(XHitSize, XHitArr, outChi2.ptr, minChi2_reg, bestHit_reg, hit_cnt, N);
  }
  // TODO
  // get the maxCandsPerSeed best cands;
  // Eventually add to the list of cands for this SEED
  // ...
}


__global__ void findCandidates_kernel(
    Hit *hits, const GPlexQI XHitSize, const GPlexHitIdx XHitArr, 
    const GPlexLS propErr, GPlexHS msErr, GPlexHV msPar,
    const GPlexLV propPar, GPlexQF outChi2,
    GPlexQF Chi2, GPlexQI HitsIdx,
    const int maxSize, const int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  findCandidates_fn(hits, XHitSize, XHitArr, 
    propErr, msErr, msPar,
    propPar, outChi2,
    Chi2, HitsIdx,
    maxSize, itrack, N);
}


void findCandidates_wrapper(const cudaStream_t &stream,
    LayerOfHitsCU &layer, const GPlexQI &XHitSize,  const GPlexHitIdx &XHitArr,
    const GPlexLS &propErr, GPlexHS &msErr, GPlexHV &msPar,
    const GPlexLV &propPar, GPlexQF &outChi2,
    GPlexQF &Chi2, GPlexQI &HitsIdx,
    const int maxSize, const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  findCandidates_kernel <<< grid, block, 0, stream >>>
    (layer.m_hits, XHitSize, XHitArr,
     propErr, msErr, msPar, propPar, outChi2,
     /*propErr.ptr, propErr.stride,*/
     /*msErr.ptr, msErr.stride, msErr.kSize,*/
     /*msPar.ptr, msPar.stride, msPar.kSize,*/
     /*outChi2.ptr, outChi2.stride,*/
     Chi2, HitsIdx,
     maxSize, N);
}
#endif
