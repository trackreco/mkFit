#include "best_hit_kernels.h"

#include "computeChi2_kernels.h"
#include "index_selection_kernels.h"
#include "reorganize_gplex.h"
#include "array_algorithms_cu.h"
#include "kalmanUpdater_kernels.h"
#include "propagation_kernels.h"

#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

#define BLOCK_SIZE_X 256


__device__ void getNewBestHitChi2_fn(
    const GPlexQI &XHitSize, const GPlexHitIdx &XHitArr,
    const float *outChi2, float &minChi2,
    int &bestHit, const int hit_cnt, const int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;

  if (itrack < N) {
    if (hit_cnt < XHitSize[itrack]) {
      float chi2 = fabs(outChi2[itrack]);
      if (chi2 < minChi2) {
        minChi2 = chi2;
        bestHit = XHitArr(itrack, hit_cnt, 0);
      }
    }
  }
}


__global__ void getNewBestHitChi2_kernel(
    const GPlexQI XHitSize, const GPlexHitIdx XHitArr,
    const float *outChi2, float *minChi2,
    int *bestHit, const int hit_cnt, const int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  if (itrack < N) {
    getNewBestHitChi2_fn(XHitSize, XHitArr, outChi2, 
        minChi2[itrack], bestHit[itrack], hit_cnt, N);
  }
}


void getNewBestHitChi2_wrapper(const cudaStream_t &stream,
    const GPlexQI &XHitSize, const GPlexHitIdx &XHitArr,
    const GPlexQF &outChi2, float *minChi2, int *bestHit, 
    const int hit_cnt, const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  getNewBestHitChi2_kernel <<< grid, block, 0, stream >>>
    (XHitSize, XHitArr, outChi2.ptr, minChi2, bestHit, hit_cnt, N);
}


void fill_array_cu(float *array, const int size, const float value) {
  thrust::device_ptr<float> d_ptr(array);
  thrust::fill(d_ptr, d_ptr + size, value);
}


__device__ void updateTracksWithBestHit_fn(Hit *hits, 
    const float minChi2, const int bestHit,
    GPlexHS &msErr, GPlexHV &msPar, const GPlexLV &propPar, 
    GPlexQF &Chi2, GPlexQI &HitsIdx, const int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  if (itrack < N) {
    if (bestHit >= 0)
    {
      Hit   &hit  = hits[ bestHit ];
      const float &chi2_local = minChi2;
	  
      for (int i = 0; i < msErr.kSize; ++i) {
        msErr[itrack + i*msErr.stride] = hit.errArrayCU()[i];
      }
      for (int i = 0; i < msPar.kSize; ++i) {
        msPar(itrack, i, 0) = hit.posArrayCU()[i];
      }
      Chi2[itrack] += chi2_local;
      HitsIdx[itrack] = bestHit;
    }
    else
    {
      msErr[itrack+ 0*msErr.stride] = 666;
      msErr[itrack+ 1*msErr.stride] = 0;
      msErr[itrack+ 2*msErr.stride] = 666;
      msErr[itrack+ 3*msErr.stride] = 0;
      msErr[itrack+ 4*msErr.stride] = 0;
      msErr[itrack+ 5*msErr.stride] = 666;

      for (int i = 0; i < msPar.kSize; ++i) {
        msPar(itrack, i, 0) = propPar(itrack, i, 0);
      }
      HitsIdx[itrack] = -1;
      // Don't update chi2
    }
  }
}


__global__ void updateTracksWithBestHit_kernel(Hit *hits, 
    const float *minChi2, const int *bestHit,
    GPlexHS msErr, GPlexHV msPar, const GPlexLV propPar, 
    GPlexQF Chi2, GPlexQI HitsIdx, const int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  if (itrack < N) {
    updateTracksWithBestHit_fn
        (hits, minChi2[itrack], bestHit[itrack],
         msErr, msPar, propPar, Chi2, HitsIdx, N);
  }
}


void updateTracksWithBestHit_wrapper(const cudaStream_t &stream,
    LayerOfHitsCU &layer, const float *minChi2, const int *best_hit, 
    GPlexHS &msErr, GPlexHV &msPar, const GPlexLV &propPar,
    GPlexQF &Chi2, GPlexQI &HitsIdx, const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  updateTracksWithBestHit_kernel <<< grid, block, 0, stream >>>
      (layer.m_hits.data(), minChi2, best_hit, msErr, msPar, propPar, Chi2, HitsIdx, N);
}


int getMaxNumHits_wrapper(const GPlexQI d_XHitSize, const int N) {
  thrust::device_ptr<int> d_ptr(d_XHitSize.ptr);
  int maxSize=  thrust::reduce(d_ptr, d_ptr + N, -1, thrust::maximum<int>());
  maxSize = std::min(maxSize, Config::maxHitsConsidered);

  return maxSize;
}


__device__ void bestHit_fn(
    Hit *hits, const GPlexQI &XHitSize, const GPlexHitIdx &XHitArr, 
    const GPlexLS &propErr, GPlexHS &msErr, GPlexHV &msPar,
    const GPlexLV &propPar, GPlexQF &outChi2,
    GPlexQF &Chi2, GPlexQI &HitsIdx,
    const int maxSize, const int itrack, const int N) {

  /*int itrack = threadIdx.x + blockDim.x*blockIdx.x;*/
  int bestHit_reg = -1;
  float minChi2_reg = Config::chi2Cut;

  if (itrack < N)
    HitsIdx[itrack] = 0;

  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
    HitToMs_fn(msErr, msPar, hits, XHitSize, XHitArr, HitsIdx, hit_cnt, itrack, N);
#if 0
      // TODO: add CMSGeom
      if (Config::useCMSGeom) {
        //propagateHelixToRMPlex(psErr,  psPar, inChg,  msPar, propErr, propPar);
        throw std::runtime_error("useCMSGeom not implemented yet for GPU");
      } else {}
#endif
    computeChi2_fn(propErr, msErr, msPar, propPar, outChi2, itrack, N);
    getNewBestHitChi2_fn(XHitSize, XHitArr, outChi2.ptr, minChi2_reg, bestHit_reg, hit_cnt, N);
  }
  updateTracksWithBestHit_fn
      (hits, 
       minChi2_reg, bestHit_reg,
       msErr, msPar, propPar,
       Chi2, HitsIdx,
       N);
}


__global__ void bestHit_kernel(
    Hit *hits, const GPlexQI XHitSize, const GPlexHitIdx XHitArr, 
    const GPlexLS propErr, GPlexHS msErr, GPlexHV msPar,
    const GPlexLV propPar, GPlexQF outChi2,
    GPlexQF Chi2, GPlexQI HitsIdx,
    const int maxSize, const int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  bestHit_fn(hits, XHitSize, XHitArr, 
    propErr, msErr, msPar,
    propPar, outChi2,
    Chi2, HitsIdx,
    maxSize, itrack, N);
}


void bestHit_wrapper(const cudaStream_t &stream,
    LayerOfHitsCU &layer, const GPlexQI &XHitSize,  const GPlexHitIdx &XHitArr,
    const GPlexLS &propErr, GPlexHS &msErr, GPlexHV &msPar,
    const GPlexLV &propPar, GPlexQF &outChi2,
    GPlexQF &Chi2, GPlexQI &HitsIdx,
    const int maxSize, const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  bestHit_kernel <<< grid, block, 0, stream >>>
    (layer.m_hits.data(), XHitSize, XHitArr,
     propErr, msErr, msPar, propPar, outChi2,
     Chi2, HitsIdx,
     maxSize, N);
}


template <int BLOCK_THREADS>
__global__ void findBestHit_kernel(LayerOfHitsCU *layers,
                                   EtaBinOfCandidatesCU *etabin_of_cands,
                                   GPlexQI XHitSize, GPlexHitIdx XHitArr,
                                   GPlexLS Err_iP, GPlexLV Par_iP,
                                   GPlexHS *msErr_arr, GPlexHV *msPar_arr,
                                   GPlexLS Err_iC, GPlexLV Par_iC,
                                   GPlexQF outChi2,
                                   GPlexQF Chi2, GPlexQI *HitsIdx_arr,
                                   GPlexQI inChg, GPlexQI Label, GeometryCU geom, 
                                   int *maxSize, int gplex_size,
                                   const int nlayers_per_seed) {
  for (int ebin = 0; ebin != Config::nEtaBin; ++ebin) {
    for (int beg = 0; beg < etabin_of_cands[ebin].m_fill_index; beg += gplex_size) {
      int end = min(beg + gplex_size, etabin_of_cands[ebin].m_fill_index);
      int N = end - beg; 

      int tidx = threadIdx.x + blockDim.x*blockIdx.x;
      int itrack = beg + tidx;

      if (itrack < end) {

        InputTracksCU_fn(etabin_of_cands[ebin].m_candidates, Err_iP, Par_iP,
            inChg, Chi2, Label, HitsIdx_arr, beg, end, tidx, N);

        for (int ilay = nlayers_per_seed; ilay < Config::nLayers; ++ilay)
        {
          int hit_idx = ilay;
          GPlexHS &msErr = msErr_arr[hit_idx];
          GPlexHV &msPar = msPar_arr[hit_idx];
          GPlexQI &HitsIdx = HitsIdx_arr[hit_idx];

          float *radii = geom.radii;

          LayerOfHitsCU &layer = layers[ilay];

          int maxSize_block;
          selectHitIndices_fn(layer, Err_iP, Par_iP, XHitSize, XHitArr, tidx, N);
          reduceMax_fn<int, BLOCK_THREADS, 1, cub::BLOCK_REDUCE_WARP_REDUCTIONS>
            (XHitSize.ptr, XHitSize.N, &maxSize_block);
          bestHit_fn(layer.m_hits.data(), XHitSize, XHitArr,
                     Err_iP, msErr, msPar, Par_iP, outChi2,
                     Chi2, HitsIdx, maxSize_block, tidx, N);
          kalmanUpdate_fn( Err_iP, msErr, Par_iP, msPar, Par_iC, Err_iC, tidx, N);
          if (ilay+1 < Config::nLayers) {
            float radius = radii[ilay+1];
            propagationForBuilding_fn(Err_iC, Par_iC, inChg, radius, Err_iP, Par_iP, tidx, N);
          }
        }
        OutputTracksCU_fn(etabin_of_cands[ebin].m_candidates, 
            Err_iP, Par_iP, inChg, Chi2, Label, HitsIdx_arr, beg, end, tidx, N);
      }
    }
  }
}


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
    int *maxSize, int N) {
  int gridx = (N-1)/BLOCK_SIZE_X + 1;
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  if (gridx > max_blocks_x) {
    throw std::runtime_error("The matriplex size should be chosen such "
                             "that gplex_size*BLOCK_SIZE_X <= max_blocks_x");
  }
  // The loop over tracks is taking care of the case where there is more tracks
  // than available threads.
  // We should actually not throw an exception, GPlex should be allocated with
  // a smaller size in MkFitter.

  findBestHit_kernel<BLOCK_SIZE_X> <<< grid, block, 0, stream >>>
    (layers, event_of_cands_cu.m_etabins_of_candidates,
     XHitSize, XHitArr, Err_iP, Par_iP, msErr, msPar,
     Err_iC, Par_iC, outChi2, Chi2, HitsIdx, inChg, Label, geom, maxSize, N,
     Config::nlayers_per_seed);
}
