#include "gpu_utils.h"

void sync_gpu() {
  cudaCheckErrorSync();
}
