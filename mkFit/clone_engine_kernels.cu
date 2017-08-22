#include "clone_engine_kernels.h"
#include "clone_engine_cuda_helpers.h"

#include "BestCands.h"

#include "computeChi2_kernels.h"
#include "index_selection_kernels.h"
#include "reorganize_gplex.h"
#include "build_tracks_kernels.h"
#include "array_algorithms_cu.h"
#include "propagation_kernels.h"


constexpr int BLOCK_SIZE_X = 256;


__device__ void findCandidates_fn(const int ilay,
    Hit *hits, const GPlexQI &XHitSize, const GPlexHitIdx &XHitArr, 
    const GPlexLS &propErr, GPlexHS &msErr, GPlexHV &msPar,
    const GPlexLV &propPar, GPlexQF &outChi2,
    GPlexQF &Chi2, GPlexQI *HitsIdx_arr,
    GPlexQB &Valid,
    const int maxSize, const int itrack, const int N,
    CandsGPU::BestCands<Config::maxCandsPerSeed, BLOCK_SIZE_X>& cand_list)
{
  if (itrack >= N) return;
  if (!Valid[itrack]) return;

  const int icand = itrack % Config::maxCandsPerSeed;

  findTrackCands_fn(ilay, hits, XHitSize, XHitArr, 
    propErr, msErr, msPar, propPar, outChi2, Chi2, HitsIdx_arr,
    maxSize, itrack, icand, N, cand_list);

  // At this point, each track has the best candidates organized in a stack
  //__syncthreads(); no need for that, sync hapeen in merge
  // sorting garbbage is fine (i.e. threads with no valid tracks to start with)
  // because of the default sentinel values.
  int itrack_block = threadIdx.x;
  int iseed_block = itrack_block / Config::maxCandsPerSeed;
  int icand_block = itrack_block % Config::maxCandsPerSeed;
  cand_list.merge_cands_for_seed(iseed_block, icand_block);
  // At this point the best overall candidates for seed iseed are
  // stored in the first list (icand == 0)
}


__device__ void updateCandsWithInfo(Track *candidates, int *ntracks_per_seed,
                                    int my_trkIdx, int my_hitIdx, int my_nhits,
                                    float my_chi2, 
                                    int beg, int itrack)
{
  int itrack_ev = beg + itrack;
  int iseed_ev = itrack_ev / Config::maxCandsPerSeed;
  int icand_ev = itrack_ev % Config::maxCandsPerSeed; 

  int ncands_seed = ntracks_per_seed[iseed_ev];

  Track tmp_track;  // BAD: too big to fit in registers

  if (0 <= my_trkIdx && my_trkIdx < ncands_seed) {
    // Get candidate of previous step that match this trkIdx
    Track& base_track = candidates[iseed_ev*Config::maxCandsPerSeed + my_trkIdx];

    tmp_track = base_track;
    tmp_track.addHitIdx(my_hitIdx, 0.f /*chi2*/);
    tmp_track.setChi2(my_chi2);
  }

  // Before writing the new track to the EtaBinOfCombCandidatesCU
  // makes sure everybody has a copy of it.
  // CUDA 8; synchronize only cands for seed;
  __syncthreads();

  if (0 <= my_trkIdx && my_trkIdx < ncands_seed) {
    Track& my_track = candidates[iseed_ev*Config::maxCandsPerSeed + icand_ev];
    my_track = tmp_track;
  }
}

 
__global__ void findInLayers_kernel(LayerOfHitsCU *layers, 
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
  constexpr int max_cands_per_block = (BLOCK_SIZE_X / Config::maxCandsPerSeed) * Config::maxCandsPerSeed;

  for (int ebin = 0; ebin < Config::nEtaBin; ebin++) {
    EtaBinOfCombCandidatesCU& etabin_of_cands = etabins_of_comb_candidates[ebin];
    int nseed = etabin_of_cands.m_nseed;
    if (nseed == 0) continue;

    int max_nseed_plex = gplex_size/Config::maxCandsPerSeed;

    for (int beg_seed = 0; beg_seed < nseed; beg_seed += max_nseed_plex) {
      int beg = beg_seed * Config::maxCandsPerSeed;
      int end_seed = min(beg_seed + max_nseed_plex, nseed);
      int end = end_seed * Config::maxCandsPerSeed;

      int itrack_plex = threadIdx.x + blockIdx.x * max_cands_per_block;  // starts at 0. Fits a gplex 
      int itrack /* itrack_ev */  = beg + itrack_plex;  // absolute track idx. If 'num_tracks' > gplex_size
      int iseed_ev = itrack / Config::maxCandsPerSeed;
      int icand_ev = itrack % Config::maxCandsPerSeed;
      int N = end - beg;

      int itrack_block = threadIdx.x;
      int iseed_block = itrack_block / Config::maxCandsPerSeed;
      int icand_block = itrack_block % Config::maxCandsPerSeed;

      int maxSize_block;

      __shared__ CandsGPU::BestCands<Config::maxCandsPerSeed, BLOCK_SIZE_X> cand_list;

      // NOTE: maxCandsPerSeed should divide BLOCK_SIZE_X because
      //       padding is handled differently in GPlex and CandList (in SM)
      //       One easy way to deal with that is to allocate GPlexes a bit wider and
      //       let the last threads of blocks being idle; and referencing arrays with:
      //            tid = tidx + bidx * BLOCK_SIZE_X;
      //       we still need itrack_plex for the guard.
      //       We might has well not do it, and take advantages of these thread to
      //       enlarge the number of candidates for each seed.
      if (itrack_block < max_cands_per_block && itrack_plex < N) {
          for (int ilay = nlayers_per_seed; ilay <= Config::nLayers; ++ilay) {
            Track *tracks = etabin_of_cands.m_candidates.data();
            int *tracks_per_seed = etabin_of_cands.m_ntracks_per_seed.data();

            int Nhits = ilay;
            InputTracksAndHitIdxComb_fn(tracks, tracks_per_seed, 
                Err_iP, Par_iP,
                inChg, Chi2, Label, HitsIdx_arr, 
                SeedIdx, CandIdx, Valid, Nhits ,
                beg, end, itrack_plex, N);

            // no cands at this position for this layer
            if (Valid[itrack_plex]) {
              if (ilay > nlayers_per_seed) {
                UpdateWithLastHit_fn(layers[ilay-1], HitsIdx_arr[ilay-1],
                    msPar_arr[ilay-1], msErr_arr[ilay-1],
                    Par_iP, Err_iP, Par_iC, Err_iC, itrack_plex, N);
                if (ilay < Config::nLayers) {
                  propagationForBuilding_fn(Err_iC, Par_iC, inChg, geom.radii[ilay], 
                      Err_iP, Par_iP, itrack_plex, N);
                  OutputParErrCU_fn(tracks, Err_iP, Par_iP,
                      beg, end, itrack_plex, N);
                } else {
                  OutputParErrCU_fn(tracks, Err_iC, Par_iC,
                      beg, end, itrack_plex, N);
                  break;
                }
              }
            } else {
              if (ilay >= Config::nLayers) {  // for those invalid hits
                break;
              }
            }
            XHitSize[itrack_plex] = 0;  // reset here has well, because of invalid tracks
            if (Valid[itrack_plex]) {
              selectHitIndices_fn(layers[ilay], Err_iP, Par_iP, XHitSize, XHitArr, itrack_plex, N);
            }
            reduceMaxPartial_fn<int, BLOCK_SIZE_X, 1, cub::BLOCK_REDUCE_WARP_REDUCTIONS>
              (XHitSize.ptr, 
               max_cands_per_block * blockIdx.x,
               max_cands_per_block,
               &maxSize_block);

            cand_list.reset(itrack_block);

            findCandidates_fn(ilay,
                layers[ilay].m_hits.data(), XHitSize, XHitArr,
                Err_iP, msErr_arr[ilay], msPar_arr[ilay], Par_iP, 
                outChi2, Chi2, HitsIdx_arr, Valid,
                maxSize_block, itrack_plex, N, cand_list);
            // cand_list: each seed has the best cands 
            // It's similar to CandCloner::ProcessSeddRange (called from end_smthg)
            // except than everything is already sorted.
            // Here we cannot avoid bank conflicts, because we scatter (as in MPI_Scatter)
            // a column (representing the list of overall best candidates) to different=
            // thread. This is a consequence of the BestCand structure that
            // was designed to minimized bank conflicts when performing the more 
            // frequent stack insertion operation.

            // Don't check if Valid as a previously inValid track might pick a correct value
            int my_seed_track_idx = iseed_block * Config::maxCandsPerSeed;
            int my_trkIdx, my_hitIdx, my_nhits;
            float my_chi2;

            cand_list.get_cand_info(my_seed_track_idx, icand_block,
                my_trkIdx, my_hitIdx, my_nhits, my_chi2);
            // Needs to count how many non sentinel there is:
            if (!icand_ev) {
              tracks_per_seed[iseed_ev] = cand_list.count_valid_cands(itrack_block);
            }
            __syncthreads();

            // Now, every thread has the info of its own new best cand.
            // -> Feed the etabin with that
            updateCandsWithInfo(tracks, tracks_per_seed, 
                my_trkIdx, my_hitIdx, my_nhits, my_chi2, 
                beg, itrack_plex);
          }  // ilay
          // final sorting, 
      }
    }
  }
}

__global__ void print_stuffs(Track* tracks, int*ntracks_per_seed) {
  if (threadIdx.x + blockIdx.x * blockDim.x == 0) {
    int idx = 0 + 1*Config::maxCandsPerSeed ;
    int val = tracks[idx].getHitIdx(0);
    printf("tracks[%d] = %d\n", idx, val);
    val = tracks[idx].getHitIdx(1);
    printf("tracks[%d] = %d\n", idx, val);
    val = tracks[idx].getHitIdx(2);
    printf("tracks[%d] = %d\n", idx, val);
    val = ntracks_per_seed[idx]; 
    printf("tracks[%d] = %d\n", idx, val);
  }
}

void findInLayers_wrapper(cudaStream_t &stream,
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

  findInLayers_kernel <<<grid, block, 0, stream>>>
      (layers, 
       event_of_cands_cu.m_etabins_of_comb_candidates.data(),
       XHitSize, XHitArr, Err_iP, Par_iP, msErr, msPar,
       Err_iC, Par_iC, outChi2, Chi2, HitsIdx_arr, inChg, Label, 
       SeedIdx, CandIdx, Valid, geom, maxSize, N, Config::nlayers_per_seed);
  
  /*cudaCheckErrorSync();*/
}
