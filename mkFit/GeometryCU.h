#ifndef GEOMETRY_CU_H
#define GEOMETRY_CU_H

struct GeometryCU {
  float *radii;

  void allocate() {
    cudaMalloc((void**)&radii, Config::nLayers * sizeof(float));
  }
  void deallocate() {
    cudaFree(radii);
  }
  void getRadiiFromCPU(const float *h_radii) {
    cudaMemcpy(radii, h_radii, Config::nLayers * sizeof(float), cudaMemcpyHostToDevice);
  }
};

#endif /* ifndef GEOMETRY_CU_H */
