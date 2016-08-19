#include "check_gpu_hit_structures.h"

/*#include "reorganize_gplex.cu"*/
#include "HitStructures.h"
#include "HitStructuresCU.h"
#include "reorganize_gplex.h"

#include <iostream>


__global__ void get_hit_pos_and_err(LayerOfHitsCU *layers,
    int ilay, int hit_idx, float *pos, float *err, int pos_size, int err_size) {
  if (threadIdx.x + blockDim.x * blockIdx.x == 0) {
    LayerOfHitsCU &layer = layers[ilay];
    Hit &hit = layer.m_hits[hit_idx];
    float *posArray = get_posArray(hit);
    float *errArray = get_errArray(hit);
    for (int i = 0; i < pos_size; ++i) {
      pos[i] = posArray[i];
    }
    for (int i = 0; i < err_size; ++i) {
      err[i] = errArray[i];
    }
  }
}


void compare_carrays(const float *h_a, const float *d_a, 
                     const float prec, const int n) 
{
  for (int i = 0; i < n; ++i) {
    // should be relative comparison, verify if div by 0 will happen
    if (std::abs(h_a[i] - d_a[i]) > prec) {
      std::cerr << i << " : " << h_a[i] << " / " << d_a[i] << std::endl;
    }
  }
}


void check_event_of_hits_gpu(const EventOfHits& event_of_hits)
{
  EventOfHitsCU event_of_hits_cu;
  event_of_hits_cu.allocGPU(event_of_hits);
  event_of_hits_cu.copyFromCPU(event_of_hits);

  constexpr int pos_size = 3;
  constexpr int err_size = 6;

  float *d_pos, *d_err;
  float pos[pos_size], err[err_size];

  cudaMalloc((void**)&d_pos, pos_size*sizeof(float));
  cudaMalloc((void**)&d_err, err_size*sizeof(float));

  dim3 grid(1, 1, 1);
  dim3 block(1, 1, 1);

  int ilay = 2;
  int hit_idx = 3;

  get_hit_pos_and_err <<< grid, block >>>
    (event_of_hits_cu.m_layers_of_hits, ilay, hit_idx, d_pos, d_err, pos_size, err_size);

  cudaMemcpy(pos, d_pos, pos_size*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(err, d_err, err_size*sizeof(float), cudaMemcpyDeviceToHost);

  compare_carrays(event_of_hits.m_layers_of_hits[ilay].m_hits[hit_idx].posArray(),
                  pos, 1e-3, pos_size);
  compare_carrays(event_of_hits.m_layers_of_hits[ilay].m_hits[hit_idx].errArray(),
                  err, 1e-3, err_size);

  cudaFree(d_pos);
  cudaFree(d_err);

  event_of_hits_cu.deallocGPU();
}
