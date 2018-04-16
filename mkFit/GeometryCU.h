#ifndef GEOMETRY_CU_H
#define GEOMETRY_CU_H

#include "gpu_utils.h"

namespace mkfit {

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

} // end namespace mkfit
#endif /* ifndef GEOMETRY_CU_H */
