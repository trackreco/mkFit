#ifndef _GPLEX_H_
#define _GPLEX_H_

#include <cuda_runtime.h>
#include <stdio.h>

template <typename T>
struct GPlex { 
  T* ptr;
  size_t pitch, stride, x, y;

  void allocate(size_t ntracks, size_t plex_size) {
    x = ntracks;
    y = plex_size;
    cudaMallocPitch((void**)&ptr, &pitch, x*sizeof(T), y);
    stride = pitch/sizeof(T);
  }
  void free() {
    cudaFree(ptr);
    x = 0; y = 0; pitch = 0; stride = 0;
  }
  //cudaMemcpy2D(d_msErr.ptr, d_msErr.pitch, msErr.fArray, N*sizeof(T),
               //N*sizeof(T), HS, cudaMemcpyHostToDevice);
};

#endif  // _GPLEX_H_
