#ifndef _GPLEX_H_
#define _GPLEX_H_

#include <cuda_runtime.h>
#include <stdio.h>

#include "gpu_utils.h"
#include "Matrix.h"

//#include "gpu_constants.h"
__device__ __constant__ static int gplexSymOffsets[7][36] =
{
  {}, 
  {},
  { 0, 1, 1, 2 },
  { 0, 1, 3, 1, 2, 4, 3, 4, 5 }, // 3
  {},
  {},
  { 0, 1, 3, 6, 10, 15, 1, 2, 4, 7, 11, 16, 3, 4, 5, 8, 12, 17, 6, 7, 8, 9, 13, 18, 10, 11, 12, 13, 14, 19, 15, 16, 17, 18, 19, 20 }
};


// GPU implementation of a Matriplex-like structure
// The number of tracks is the fast dimension and is padded in order to have
// consecutive and aligned memory accesses. For cached reads, this result in a
// single memory transaction for the 32 threads of a warp to access 32 floats.
// See:
// http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#global-memory-3-0
// In practice, The number of tracks (ntracks) is set to be MPT_SIZE
template <typename M>
struct GPlex { 
  using T = typename M::value_type;
  using value_type = T;

  size_t kRows = M::kRows;
  size_t kCols = M::kCols;
  size_t kSize = M::kSize;

  T* ptr;
  size_t pitch, stride, N;

  __device__ T  operator[](int xx) const { return ptr[xx]; }
  __device__ T& operator[](int xx)       { return ptr[xx]; }

  __device__ T& operator()(int n, int i, int j)       { return ptr[n + (i*kCols + j)*stride]; }
  __device__ T  operator()(int n, int i, int j) const { return ptr[n + (i*kCols + j)*stride]; }

  void allocate(size_t ntracks) {
    N = ntracks;
    cudaMallocPitch((void**)&ptr, &pitch, N*sizeof(T), kSize);
    stride = pitch/sizeof(T);  // Number of elements
  }
  void free() {
    cudaFree(ptr);
    N = 0; pitch = 0; stride = 0;
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
  void copyAsyncFromDevice(cudaStream_t& stream, GPlex<M>& gplex) {
    cudaMemcpy2DAsync(ptr, pitch, gplex.ptr, gplex.pitch,
                      N*sizeof(T), kSize, cudaMemcpyDeviceToDevice, stream);
    cudaCheckError();
  }
};


template <typename M>
struct GPlexSym : GPlex<M> {
  using T = typename GPlex<M>::T;
  using GPlex<M>::kRows;
  using GPlex<M>::kCols;
  using GPlex<M>::stride;
  using GPlex<M>::ptr;
  __device__ size_t Off(size_t i) const { return gplexSymOffsets[kRows][i]; }
  // Note: convenient but noticeably slower due to the indirection
  __device__ T& operator()(int n, int i, int j)       { return ptr[n + Off(i*kCols + j)*stride]; }
  __device__ T  operator()(int n, int i, int j) const { return ptr[n + Off(i*kCols + j)*stride]; }
  //__device__ T& operator()(int n, int i, int j)       { return ptr[n + i*stride]; }
  //__device__ T  operator()(int n, int i, int j) const { return ptr[n + i*stride]; }
};

using GPlexLL = GPlex<MPlexLL>;
using GPlexLV = GPlex<MPlexLV>;
using GPlexLS = GPlexSym<MPlexLS>;

using GPlexHH = GPlex<MPlexHH>;
using GPlexHV = GPlex<MPlexHV>;
using GPlexHS = GPlexSym<MPlexHS>;

using GPlexLH = GPlex<MPlexLH>;

using GPlexQF = GPlex<MPlexQF>;
using GPlexQI = GPlex<MPlexQI>;
using GPlexQB = GPlex<MPlexQB>;

const int GPlexHitIdxMax = 16;  // FIXME: copied and past from MkFitter.h
using GPlexHitIdx = GPlex<Matriplex::Matriplex<int, GPlexHitIdxMax, 1, MPT_SIZE>>;

template <typename M>
struct GPlexReg {
  using T = typename M::value_type;
  using value_type = T;

  size_t kRows = M::kRows;
  size_t kCols = M::kCols;
  size_t kSize = M::kSize;

  __device__ T  operator[](int xx) const { return arr[xx]; }
  __device__ T& operator[](int xx)       { return arr[xx]; }

  __device__ T& operator()(int n, int i, int j)       { return arr[i*kCols + j]; }
  __device__ T  operator()(int n, int i, int j) const { return arr[i*kCols + j]; }

  __device__ void SetVal(T v)
  {
     for (int i = 0; i < kSize; ++i)
     {
        arr[i] = v;
     }
  }

  T arr[M::kSize];
};

using GPlexRegLL = GPlexReg<MPlexLL>;
using GPlexRegLH = GPlexReg<MPlexLH>;
using GPlexRegHH = GPlexReg<MPlexHH>;
using GPlexRegLV = GPlexReg<MPlexLV>;
using GPlexRegHS = GPlexReg<MPlexHS>;
using GPlexRegHV = GPlexReg<MPlexHV>;
using GPlexReg2V = GPlexReg<MPlex2V>;
using GPlexReg2S = GPlexReg<MPlex2S>;
using GPlexRegQF = GPlexReg<MPlexQF>;

#endif  // _GPLEX_H_
