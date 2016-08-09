#ifndef GPU_UTILS_H
#define GPU_UTILS_H 

#include <cub/util_debug.cuh>

#define cudaCheckError()               \
  do {                                 \
    cudaError_t e=cudaGetLastError();  \
    CubDebugExit(e)                    \
  } while(0)  

#define cudaCheckErrorSync()           \
  do {                                 \
    cudaDeviceSynchronize();           \
    cudaCheckError();                  \
  } while(0)

// CUDA specific:
// Maximum number of blocks in the X direction of the thread grid.
constexpr int max_blocks_x = 1 << 15;

// The first call to a CUDA API function takes the initialization hit.
void separate_first_call_for_meaningful_profiling_numbers();

void sync_gpu();

#endif /* ifndef GPU_UTILS_H */
