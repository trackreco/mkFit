#ifndef GEOMETRY_CU_H
#define GEOMETRY_CU_H

#include "gpu_utils.h"

struct GeometryCU {
  float *radii = nullptr;

  void allocate() {
    if (radii == nullptr)
      cudaMalloc((void**)&radii, Config::nLayers * sizeof(float));
    //cudaCheckError();
  }
  void deallocate() {
    if (radii != nullptr) {
      cudaFree(radii);
      cudaCheckError();
      radii = nullptr;
    }
  }
  void getRadiiFromCPU(const float *h_radii) {
    cudaMemcpy(radii, h_radii, Config::nLayers * sizeof(float), cudaMemcpyHostToDevice);
    //cudaCheckError();
  }
};

#endif /* ifndef GEOMETRY_CU_H */
