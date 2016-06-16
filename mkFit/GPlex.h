#ifndef _GPLEX_H_
#define _GPLEX_H_

#include <cuda_runtime.h>
#include <stdio.h>

#include "Matrix.h"

#define cudaCheckError() {                                          \
  cudaError_t e=cudaGetLastError();                                 \
  if(e!=cudaSuccess) {                                              \
    printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
    exit(0); \
  }                                                                 \
}

// GPU implementation of a Matriplex-like structure
// The number of tracks is the fast dimension and is padded in order to have
// consecutive and aligned memory accesses. For cached reads, this result in a
// single memory transaction for the 32 threads of a warp to access 32 floats.
// See:
// http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#global-memory-3-0
// In practice, The number of tracks (ntracks) is set to be MPT_SIZE
template <typename T, typename M>
struct GPlex { 
  using value_type = T;
  T* ptr;
  size_t pitch, stride, N, kSize;

  __device__ T  operator[](int xx) const { return ptr[xx]; }
  __device__ T& operator[](int xx)       { return ptr[xx]; }

  __device__ T& operator()(int n, int i, int j)       { return ptr[n + i*stride]; }
  __device__ T  operator()(int n, int i, int j) const { return ptr[n + i*stride]; }

  void allocate(size_t ntracks, size_t aSize) {
    N = ntracks;
    kSize = aSize;
    cudaMallocPitch((void**)&ptr, &pitch, N*sizeof(T), kSize);
    stride = pitch/sizeof(T);  // Number of elements
  }
  void free() {
    cudaFree(ptr);
    N = 0; kSize = 0; pitch = 0; stride = 0;
  }
  //cudaMemcpy2D(d_msErr.ptr, d_msErr.pitch, msErr.fArray, N*sizeof(T),
               //N*sizeof(T), HS, cudaMemcpyHostToDevice);

  void copyAsyncFromHost(cudaStream_t& stream, const M& mplex) {
    cudaMemcpy2DAsync(ptr, pitch, mplex.fArray, N*sizeof(T),
                      N*sizeof(T), kSize, cudaMemcpyHostToDevice, stream);
    cudaCheckError();
  }
  void copyAsyncToHost(cudaStream_t& stream, M& mplex) {
    cudaMemcpy2DAsync(mplex.fArray, N*sizeof(T), ptr, pitch,
                      N*sizeof(T), kSize, cudaMemcpyDeviceToHost, stream);
    cudaCheckError();
  }
  void copyAsyncFromDevice(cudaStream_t& stream, GPlex<T, M>& gplex) {
    cudaMemcpy2DAsync(ptr, pitch, gplex.ptr, gplex.pitch,
                      N*sizeof(T), kSize, cudaMemcpyDeviceToDevice, stream);
    cudaCheckError();
  }
};

using GPlexLL = GPlex<float, MPlexLL>;
using GPlexLV = GPlex<float, MPlexLV>;
using GPlexLS = GPlex<float, MPlexLS>;

using GPlexHH = GPlex<float, MPlexHH>;
using GPlexHV = GPlex<float, MPlexHV>;
using GPlexHS = GPlex<float, MPlexHS>;

using GPlexLH = GPlex<float, MPlexLH>;

using GPlexQF = GPlex<float, MPlexQF>;
using GPlexQI = GPlex<int, MPlexQI>;

const int GPlexHitIdxMax = 16;  // FIXME: copied and past from MkFitter.h
using GPlexHitIdx = GPlex<int, Matriplex::Matriplex<int, GPlexHitIdxMax, 1, MPT_SIZE>>;

template <typename T, int D1, int D2>
struct GPlexReg {
  __device__ T  operator[](int xx) const { return arr[xx]; }
  __device__ T& operator[](int xx)       { return arr[xx]; }

  __device__ T& operator()(int n, int i, int j)       { return arr[i*D2 + j]; }
  __device__ T  operator()(int n, int i, int j) const { return arr[i*D2 + j]; }

  __device__ void SetVal(T v)
  {
     for (int i = 0; i < D1; ++i)
     {
        arr[i] = v;
     }
  }

  T arr[D1];
};

using GPlexRegLL = GPlexReg<float, 36, 6>;
using GPlexRegLH = GPlexReg<float, 18, 6>;
using GPlexRegHH = GPlexReg<float, 9, 3>;
using GPlexRegLV = GPlexReg<float, 6, 1>;
using GPlexRegHS = GPlexReg<float, 6, 1>;
using GPlexRegHV = GPlexReg<float, 3, 1>;
using GPlexReg2V = GPlexReg<float, 2, 1>;
using GPlexReg2S = GPlexReg<float, 3, 1>;
using GPlexRegQF = GPlexReg<float, 1, 1>;

#endif  // _GPLEX_H_
