#ifndef _GPLEX_H_
#define _GPLEX_H_

#include <cuda_runtime.h>
#include <stdio.h>

// GPU implementation of a Matriplex-like structure
// The number of tracks is the fast dimension and is padded in order to have
// consecutive and aligned memory accesses. For cached reads, this result in a
// single memory transaction for the 32 threads of a warp to access 32 floats.
// See:
// http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#global-memory-3-0
// In practice, The number of tracks (ntracks) is set to be MPT_SIZE
template <typename T>
struct GPlex { 
  T* ptr;
  size_t pitch, stride, x, y;

  void allocate(size_t ntracks, size_t plex_size) {
    x = ntracks;
    y = plex_size;
    cudaMallocPitch((void**)&ptr, &pitch, x*sizeof(T), y);
    stride = pitch/sizeof(T);  // Number of elements
  }
  void free() {
    cudaFree(ptr);
    x = 0; y = 0; pitch = 0; stride = 0;
  }
  //cudaMemcpy2D(d_msErr.ptr, d_msErr.pitch, msErr.fArray, N*sizeof(T),
               //N*sizeof(T), HS, cudaMemcpyHostToDevice);
};

#endif  // _GPLEX_H_
