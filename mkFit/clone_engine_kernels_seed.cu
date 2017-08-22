#include "clone_engine_kernels_seed.h"
#include "clone_engine_cuda_helpers.h"

#include "BestCands.h"

#include "computeChi2_kernels.h"
#include "index_selection_kernels.h"
#include "reorganize_gplex.h"
#include "build_tracks_kernels.h"
#include "array_algorithms_cu.h"
#include "propagation_kernels.h"

#include <cstdio>

constexpr int BLOCK_SIZE_X = 256;


__device__ void findCandidates_fn_seed(const int ilay,
    Hit *hits, const GPlexQI &XHitSize, const GPlexHitIdx &XHitArr, 
    const GPlexLS &propErr, GPlexHS &msErr, GPlexHV &msPar,
    const GPlexLV &propPar, GPlexQF &outChi2,
    GPlexQF &Chi2, GPlexQI *HitsIdx_arr,
    GPlexQB &Valid,
    const int maxSize, const int itrack, const int icand, const int N,
    CandsGPU::BestCands<Config::maxCandsPerSeed, BLOCK_SIZE_X>& cand_list)
{
  if (itrack >= N) return;

  findTrackCands_fn(ilay, hits, XHitSize, XHitArr, 
    propErr, msErr, msPar, propPar, outChi2, Chi2, HitsIdx_arr,
    maxSize, itrack, icand, N, cand_list);
}


__device__ void updateCandsWithInfo_seed(Track *candidates, int *ntracks_per_seed,
    int iseed_ev, CandsGPU::BestCands<Config::maxCandsPerSeed, BLOCK_SIZE_X>& cand_list)
{
  int iseed_block = threadIdx.x;
  Track new_cands[Config::maxCandsPerSeed];  // too big too be __shared__

  for (int icand = 0; icand < ntracks_per_seed[iseed_ev]; ++icand) {
    Track &tmp_track = new_cands[icand];

    int my_trkIdx, my_hitIdx, my_nhits;
    float my_chi2;

    cand_list.get_cand_info(iseed_block, icand, my_trkIdx, 
                            my_hitIdx, my_nhits, my_chi2);

    if (0 <= my_trkIdx && my_trkIdx < ntracks_per_seed[iseed_ev]) {
      Track& base_track = candidates[iseed_ev*Config::maxCandsPerSeed + my_trkIdx];

      tmp_track = base_track;
      tmp_track.addHitIdx(my_hitIdx, 0.f /*chi2*/);
      tmp_track.setChi2(my_chi2);
    }
  }  // icand

  // Second loop required: several new candidates can come from the same old one.
  for (int icand = 0; icand < ntracks_per_seed[iseed_ev]; ++icand) {
    candidates[iseed_ev*Config::maxCandsPerSeed + icand] = new_cands[icand];
  }
}

 
__global__ void findInLayers_kernel_seed (LayerOfHitsCU *layers, 
    EtaBinOfCombCandidatesCU* etabins_of_comb_candidates, 
    GPlexQI XHitSize, GPlexHitIdx XHitArr,  
    GPlexLS Err_iP, GPlexLV Par_iP, 
    GPlexHS *msErr_arr, GPlexHV *msPar_arr,
    GPlexLS Err_iC, GPlexLV Par_iC,
    GPlexQF outChi2,
    GPlexQF Chi2, GPlexQI *HitsIdx_arr,
    GPlexQI inChg, GPlexQI Label, 
    GPlexQI SeedIdx, GPlexQI CandIdx,
    GPlexQB Valid, GeometryCU geom, 
    int *maxSize, int gplex_size,
    const int nlayers_per_seed)
{
  int iseed_plex = threadIdx.x + blockIdx.x * blockDim.x;

  for (int ebin = 0; ebin < Config::nEtaBin; ebin++) {
    EtaBinOfCombCandidatesCU& etabin_of_cands = etabins_of_comb_candidates[ebin];
    int nseed = etabin_of_cands.m_nseed;
    if (nseed == 0) continue;

    //int max_nseed_plex = gplex_size; // no more /Config::maxCandsPerSeed; -> one thread 
    int N = gplex_size;
    int maxSize_block;

    // we need both iseed_plex and iseed_ev:
    //   * iseed_plex to access the matriplex structures, used for the algorithm
    //   * iseed_ev to access the other structures (tracks, mostly), used to 
    //     interact with the outside world.
    for (int iseed_ev = iseed_plex; 
             iseed_ev < nseed; 
             iseed_ev += gridDim.x * blockDim.x) {
      int iseed_block = threadIdx.x;

      __shared__ CandsGPU::BestCands<Config::maxCandsPerSeed, BLOCK_SIZE_X> cand_list;

      Track *tracks = etabin_of_cands.m_candidates.data();
      int *tracks_per_seed = etabin_of_cands.m_ntracks_per_seed.data();

      for (int ilay = nlayers_per_seed; ilay <= Config::nLayers; ++ilay) {
        cand_list.reset(iseed_block);
        int Nhits = ilay;

        for (int icand_seed = 0; icand_seed < tracks_per_seed[iseed_ev]; ++icand_seed) {
          InputTracksAndHitIdxComb_fn_seed(tracks, tracks_per_seed, 
              Err_iP, Par_iP,
              inChg, Chi2, Label, HitsIdx_arr, 
              SeedIdx, CandIdx, Valid, Nhits ,
              iseed_ev, icand_seed, iseed_plex, N);

          // no cands at this position for this layer
          if (ilay > nlayers_per_seed) {
            UpdateWithLastHit_fn(layers[ilay-1], HitsIdx_arr[ilay-1],
                msPar_arr[ilay-1], msErr_arr[ilay-1],
                Par_iP, Err_iP, Par_iC, Err_iC, iseed_plex, N);
            if (ilay < Config::nLayers) {
              propagationForBuilding_fn(Err_iC, Par_iC, inChg, geom.radii[ilay], 
                  Err_iP, Par_iP, iseed_plex, N);
              OutputParErrCU_fn_seed(tracks, Err_iP, Par_iP,
                  iseed_ev, icand_seed, N);
            } else {
              OutputParErrCU_fn_seed(tracks, Err_iC, Par_iC,
                  iseed_ev, icand_seed, N);
              return;
            }
          }
          XHitSize[iseed_plex] = 0;
          selectHitIndices_fn(layers[ilay], Err_iP, Par_iP, 
                              XHitSize, XHitArr, iseed_plex, N);
          reduceMaxPartial_fn<int, BLOCK_SIZE_X, 1, cub::BLOCK_REDUCE_WARP_REDUCTIONS>
            (XHitSize.ptr, blockDim.x * blockIdx.x, blockDim.x, &maxSize_block);

          findCandidates_fn_seed(ilay,
              layers[ilay].m_hits.data(), XHitSize, XHitArr,
              Err_iP, msErr_arr[ilay], msPar_arr[ilay], Par_iP, 
              outChi2, Chi2, HitsIdx_arr, Valid,
              maxSize_block, iseed_plex, icand_seed, N, cand_list);
        }  // icand_seed 
        
        cand_list.heap_sort(iseed_block, Config::maxCandsPerSeed);
        tracks_per_seed[iseed_ev] = cand_list.count_valid_cands(iseed_block);

        updateCandsWithInfo_seed(tracks, tracks_per_seed, iseed_ev, cand_list);
      }  // ilay
    }  // iseed
  }  // eta
}


void findInLayers_wrapper_seed(cudaStream_t &stream,
    LayerOfHitsCU *layers, 
    EventOfCombCandidatesCU &event_of_cands_cu,
    GPlexQI &XHitSize, GPlexHitIdx &XHitArr,
    GPlexLS &Err_iP, GPlexLV &Par_iP,
    GPlexHS *msErr, GPlexHV *msPar,
    GPlexLS &Err_iC, GPlexLV &Par_iC,
    GPlexQF &outChi2,
    GPlexQF &Chi2, GPlexQI *HitsIdx_arr,
    GPlexQI &inChg, GPlexQI &Label,
    GPlexQI &SeedIdx, GPlexQI &CandIdx,
    GPlexQB& Valid,
    GeometryCU &geom, 
    int *maxSize, int N)
{
  int gridx = (N-1)/BLOCK_SIZE_X + 1;
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  findInLayers_kernel_seed <<<grid, block, 0, stream>>>
      (layers, 
       event_of_cands_cu.m_etabins_of_comb_candidates.data(),
       XHitSize, XHitArr, Err_iP, Par_iP, msErr, msPar,
       Err_iC, Par_iC, outChi2, Chi2, HitsIdx_arr, inChg, Label, 
       SeedIdx, CandIdx, Valid, geom, maxSize, N, Config::nlayers_per_seed);
  
  /*cudaCheckErrorSync();*/
}
