#include "reorganize.h"
#include <stdio.h>

__global__ void toMatriplex_kernel(float *dst, int dst_stride,
                                   const float* __restrict__ src, int src_stride,
                                   int N, int LS) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;

  if (i < N && j < LS) {
    if (i==-1) {
      printf(" %d, mplex[%f]  /  lin[%f]\n", j, dst[i+j*dst_stride], src[j+i*src_stride]);
    }
    dst[i + j*dst_stride] = src[j + i*src_stride];
    /*float diff = fabs((dst[i + j*dst_stride] - src[j + i*src_stride]));*/
    /*if (diff > 1e-3) printf("%f\n", diff);*/
  }
}

void toMatriplex_wrapper(cudaStream_t& stream, GPlex<float> &dst, GPlex<float> &src, int N, int LS) {
  dim3 block(16, 8, 1);
  dim3 grid((N-1)/16 + 1, (LS-1)/8 +1, 1);
  toMatriplex_kernel <<<grid, block, 0, stream>>> (dst.ptr, dst.stride, src.ptr, src.stride, N, LS);
}


__global__ void reorganizeMs(float *msPar, size_t msPar_stride,
                             float *full_posArray,
                             float *msErr, size_t msErr_stride,
                             float *full_errArray,
                             int *full_hitIdx, int hi,
                             int maxHits,
                             int N, int HS, int HV, int Nhits) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;

  if (i < N) {
    int hidx = full_hitIdx[i + hi*N];
    if (j < HV) {
      /*float tmp1 = msPar[i + msPar_stride*j];*/
      msPar[i + msPar_stride*j] = full_posArray[j + HV*(hidx + hi*maxHits)];
      /*float tmp2 = msPar[i + msPar_stride*j];*/
      
      /*if (i==0 && hi == 0) {*/
        /*if (fabs(tmp1 - tmp2) > 1e-3) {*/
          /*printf("i %d, j %d, old: %f, new %f\n", i, j, tmp1, tmp2);*/
        /*}*/
      /*}*/
    }
    if (j < HS) {
      msErr[i + msErr_stride*j] = full_errArray[j + HS*(hidx + hi*maxHits)];
    }
  }
}

void reorganizeMs_wrapper(cudaStream_t& stream, GPlex<float>& msPar, float *full_posArray,
    GPlex<float>& msErr, float *full_errArray, int *full_hitIdx, int hi, int maxHits,
    int N, int hs, int hv, int Nhits) {
  dim3 block(16, 6, 1);
  dim3 grid((N-1)/16 + 1, (hs-1)/6 +1, 1);
  reorganizeMs <<<grid, block, 0, stream>>> (msPar.ptr, msPar.stride, full_posArray,
      msErr.ptr, msErr.stride, full_errArray, full_hitIdx, hi, maxHits, N, hs, hv, Nhits);
}
