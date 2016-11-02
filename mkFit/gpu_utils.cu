#include "gpu_utils.h"

void sync_gpu() {
  cudaCheckErrorSync();
}

void separate_first_call_for_meaningful_profiling_numbers() {
  sync_gpu();
}
