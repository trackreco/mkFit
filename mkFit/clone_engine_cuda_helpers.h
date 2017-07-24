#ifndef CLONE_ENGINE_CUDA_HELPERS_H
#define CLONE_ENGINE_CUDA_HELPERS_H

#include "GPlex.h"
#include "BestCands.h"
#include "reorganize_gplex.h"
#include "computeChi2_kernels.h"

__device__ int countValidHits_fn(const int itrack, const int end_hit,
                                 const GPlexQI *HitsIdx_arr);

__device__ int countInvalidHits_fn(const int itrack, const int end_hit, 
                                   const GPlexQI *HitsIdx_arr);

__device__ int2 countHits_fn(const int itrack, const int end_hit, 
                            const GPlexQI *HitsIdx_arr);

template <int MAX_CANDS_PER_SEED, int BLOCK_SIZE>
__device__ void findTrackCands_fn(const int ilay,
    Hit *hits, const GPlexQI &XHitSize, const GPlexHitIdx &XHitArr, 
    const GPlexLS &propErr, GPlexHS &msErr, GPlexHV &msPar,
    const GPlexLV &propPar, GPlexQF &outChi2,
    GPlexQF &Chi2, GPlexQI *HitsIdx_arr,
    const int maxSize, const int itrack, const int icand, const int N,
    CandsGPU::BestCands<MAX_CANDS_PER_SEED, BLOCK_SIZE>& cand_list)
{
  // Note: cand_list should initially be filled with sentinels
  int itrack_block = threadIdx.x;  // for SM arrays
  float chi2_track = Chi2[itrack];
  int XHitSize_track = XHitSize[itrack];

  GPlexQI &HitsIdx = HitsIdx_arr[ilay];

  // Note: 0 < itrack < N = matriplex_width
  HitsIdx[itrack] = 0;  // reset

  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt) {
    HitToMs_fn(msErr, msPar, hits, XHitSize, XHitArr, 
               HitsIdx, hit_cnt, itrack, N);
    // TODO: add CMSGeom
    computeChi2_fn(propErr, msErr, msPar, propPar, outChi2, itrack, N);

    // make sure the hit was in the compatiblity window for the candidate
    if (hit_cnt < XHitSize_track)
    {
      float chi2 = fabs(outChi2[itrack]); //fixme negative chi2 sometimes...
      if (chi2 < Config::chi2Cut)
      {
        int cand_hitIdx = XHitArr(itrack, hit_cnt, 0);
        int cand_nhits = countValidHits_fn(itrack, ilay, HitsIdx_arr) + 1; 
        float cand_chi2 = chi2_track + chi2;
        cand_list.update(itrack_block, icand, cand_hitIdx,
                         cand_nhits, cand_chi2);
      }
    }
  }
  // now add invalid hit 
  // if icand has enough 'good' hits, update will discard invalid hit.
  auto num_good_bad_hits = countHits_fn(itrack, ilay, HitsIdx_arr);
  int hit_idx = num_good_bad_hits.y < Config::maxHolesPerCand ? -1 : -2;
  cand_list.update(itrack_block, icand, hit_idx, 
                   num_good_bad_hits.x, 
                   chi2_track);
}

#endif  // CLONE_ENGINE_CUDA_HELPERS_H
