#include "fittracks_kernels.h"

#include "kalmanUpdater_kernels.h"
#include "propagation_kernels.h"

constexpr int BLOCK_SIZE_X = 256;

__global__ void fittracks_kernel(
      GPlexLV par_iP, GPlexLS Err_iP,
      GPlexHV *msPar_arr, GPlexHS *msErr_arr,
      GPlexLV par_iC, GPlexLS Err_iC,
      GPlexLL errorProp, GPlexQI inChg,
      const bool useParamBfield,
      const int Nhits, int N)
{
  int grid_width = blockDim.x * gridDim.x;
  int n = threadIdx.x + blockIdx.x * blockDim.x;
  for (int hi = 0; hi < Nhits; ++hi) {
    GPlexHV &msPar = msPar_arr[hi];
    GPlexHS &msErr = msErr_arr[hi];

    for (int z = 0; z < (N-1)/grid_width  +1; z++) {
      n += z*grid_width;

      propagation_fn(Err_iC, par_iC, inChg, msPar, Err_iP, par_iP, useParamBfield, n, N);
      kalmanUpdate_fn(Err_iP, msErr, par_iP, msPar, par_iC, Err_iC, n, N);
    }
  }
}

void fittracks_wrapper(cudaStream_t &stream,
                       GPlexLS &Err_iP, GPlexLV &par_iP, 
                       GPlexHS *msErr_arr, GPlexHV *msPar_arr,
                       GPlexLS &Err_iC, GPlexLV &par_iC,
                       GPlexLL &errorProp, GPlexQI &inChg,
                       const bool useParamBfield,
                       const int Nhits, const int N)
{
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  fittracks_kernel <<< grid, block, 0, stream >>>
      (par_iP, Err_iP, 
       msPar_arr, msErr_arr,
       par_iC, Err_iC,
       errorProp, inChg, 
       useParamBfield,
       Nhits, N);
}
