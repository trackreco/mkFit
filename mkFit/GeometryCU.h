#ifndef GEOMETRY_CU_H
#define GEOMETRY_CU_H

#include "gpu_utils.h"

struct GeometryCU {
  float *radii;

  void allocate() {
    cudaMalloc((void**)&radii, Config::nLayers * sizeof(float));
    cudaCheckError();
  }
  void deallocate() {
    cudaFree(radii);
    cudaCheckError();
  }
  void getRadiiFromCPU(const float *h_radii) {
    cudaMemcpy(radii, h_radii, Config::nLayers * sizeof(float), cudaMemcpyHostToDevice);
    cudaCheckError();
  }
};

#endif /* ifndef GEOMETRY_CU_H */
