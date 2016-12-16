#ifndef ARRAY_ALGORITHMS_CU_H
#define ARRAY_ALGORITHMS_CU_H 

#include <cub/cub.cuh>
//#include "atomic_utils.h"
#include "gpu_utils.h"

template <typename T,
          int BLOCK_THREADS,
          int ITEMS_PER_THREAD,
          cub::BlockReduceAlgorithm ALGORITHM>
__device__ void reduceMax_fn(const T *d_in, const int in_size, T *d_max) {
  // Specialize BlockReduce type for our thread block
  typedef cub::BlockReduce<T, BLOCK_THREADS, ALGORITHM> BlockReduceT;
  // Shared memory
  __shared__ typename BlockReduceT::TempStorage temp_storage;

  for (int i = threadIdx.x + blockIdx.x*blockDim.x;
           i < in_size;
           i += blockDim.x * gridDim.x) {
    // Per-thread tile data
    T data[ITEMS_PER_THREAD];
    int block_offset = i - threadIdx.x;
    cub::LoadDirectStriped<BLOCK_THREADS>(threadIdx.x, 
                                          d_in + block_offset, //blockIdx.x*blockDim.x,
                                          data);
    // Compute sum for a single thread block
    T aggregate = BlockReduceT(temp_storage).Reduce(data, cub::Max());

    // FIXME: Is reduction over block enough, or do we need device-wise reduction
    //        CPU code reduces on NN (typically 8 or 16) values. So block-wise should
    //        be good enough.
    //        Device-wise reductions with atomics are performance killers
    //if (threadIdx.x == 0) {
      //AtomicMax(d_max, aggregate);
    //}
    *d_max = aggregate;
  }
}

template <typename T,
          int BLOCK_THREADS,
          int ITEMS_PER_THREAD,
          cub::BlockReduceAlgorithm ALGORITHM>
__global__ void reduceMax_kernel(T *d_in, int in_size, T *d_max) {
  reduceMax_fn<T, BLOCK_THREADS, ITEMS_PER_THREAD, ALGORITHM>(d_in, in_size, d_max);
}

template <typename T, 
          int BLOCK_THREADS,
          int ITEMS_PER_THREAD>
void max_element_wrapper(T *d_v, int num_elems, T *d_max) {
  dim3 block (BLOCK_THREADS, 1, 1);
  dim3 grid (std::min(max_blocks_x, (num_elems-1)/BLOCK_THREADS + 1), 1, 1) ;
  reduceMax_kernel<T, BLOCK_THREADS, ITEMS_PER_THREAD, cub::BLOCK_REDUCE_RAKING>
       <<< grid, block >>>
       (d_v, num_elems, d_max);
}

template <typename T,
          int BLOCK_THREADS,
          int ITEMS_PER_THREAD,
          cub::BlockReduceAlgorithm ALGORITHM>
__device__ void reduceMaxPartial_fn(const T *d_in, const int start_idx,
    const int in_size, T *d_max) {
  // Specialize BlockReduce type for our thread block
  typedef cub::BlockReduce<T, BLOCK_THREADS, ALGORITHM> BlockReduceT;
  // Shared memory
  __shared__ typename BlockReduceT::TempStorage temp_storage;

  int i = threadIdx.x;
  if (i < in_size) {
    // Per-thread tile data
    T data[ITEMS_PER_THREAD];
    cub::LoadDirectStriped<BLOCK_THREADS>(threadIdx.x, 
                                          d_in + start_idx,
                                          data);
    // Compute sum for a single thread block
    T aggregate = BlockReduceT(temp_storage).Reduce(data, cub::Max());
    *d_max = aggregate;
  }
}
#endif /* ifndef ARRAY_ALGORITHMS_CU_H */
