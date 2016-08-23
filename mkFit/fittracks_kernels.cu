#include "fittracks_kernels.h"

#include "kalmanUpdater_kernels.h"
#include "propagation_kernels.h"

constexpr int BLOCK_SIZE_X = 256;

__global__ void fittracks_kernel(
      GPlexLV par_iP, GPlexLS Err_iP,
      GPlexHV msPar, GPlexHS msErr,
      GPlexLV par_iC, GPlexLS Err_iC,
      GPlexLL errorProp, GPlexQI inChg,
      int N)
{
  int grid_width = blockDim.x * gridDim.x;
  int n = threadIdx.x + blockIdx.x * blockDim.x;
  for (int z = 0; z < (N-1)/grid_width  +1; z++) {
    n += z*grid_width;

    propagation_fn(msPar, par_iC, inChg, par_iP, errorProp, Err_iP, n, N);
    kalmanUpdate_fn(Err_iP, msErr, par_iP, msPar, par_iC, Err_iC, n, N);
  }
}

void fittracks_wrapper(cudaStream_t &stream,
                       GPlexLS &Err_iP, GPlexLV &par_iP, 
                       GPlexHS *msErr, GPlexHV *msPar,
                       GPlexLS &Err_iC, GPlexLV &par_iC,
                       GPlexLL &errorProp, GPlexQI &inChg,
                       const int hit_idx, const int N)
{
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  fittracks_kernel <<< grid, block, 0, stream >>>
      (par_iP, Err_iP, 
       msPar[hit_idx], msErr[hit_idx],
       par_iC, Err_iC,
       errorProp, inChg, 
       N);
  /*kalmanUpdate_wrapper(stream, Err_iP, msErr[hit_idx],*/
                       /*par_iP, msPar[hit_idx], par_iC, Err_iC, N);*/
}
