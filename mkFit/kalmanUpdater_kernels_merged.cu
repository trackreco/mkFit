#include "Config.h"
#include "kalmanUpdater_kernels.h"

// TODO: Clean all the hard-coded #define
#define LS 21
#define HS 6
#define LH 18

#define BLOCK_SIZE_X 32
#define BLOCK_SIZE_X_UP_PARAM 32
#define BLOCK_SIZE_X_MULT_RES_ADD 256
#define MAX_BLOCKS_X 65535 // CUDA constraint


/// MultKalmanGain ////////////////////////////////////////////////////////////

__device__ void upParam_MultKalmanGain_raw(
    const float* __restrict__ a, size_t aN,
    float* sh_b, size_t bN, float *c, size_t cN, int N) {

  int n = threadIdx.x + blockIdx.x * blockDim.x;
  float sh_a[3];

  int j;
  j = 0; sh_a[0] = a[n + j*aN];
  j = 1; sh_a[1] = a[n + j*aN];
  j = 3; sh_a[2] = a[n + j*aN];

  c[ 0*cN+n] = sh_a[ 0]*sh_b[ 0] + sh_a[ 1]*sh_b[ 1] + sh_a[ 2]*sh_b[ 3];
  c[ 1*cN+n] = sh_a[ 0]*sh_b[ 1] + sh_a[ 1]*sh_b[ 2] + sh_a[ 2]*sh_b[ 4];
  c[ 2*cN+n] = sh_a[ 0]*sh_b[ 3] + sh_a[ 1]*sh_b[ 4] + sh_a[ 2]*sh_b[ 5];

  j = 1; sh_a[0] = a[n + j*aN];
  j = 2; sh_a[1] = a[n + j*aN];
  j = 4; sh_a[2] = a[n + j*aN];

  c[ 3*cN+n] = sh_a[ 0]*sh_b[ 0] + sh_a[ 1]*sh_b[ 1] + sh_a[ 2]*sh_b[ 3];
  c[ 4*cN+n] = sh_a[ 0]*sh_b[ 1] + sh_a[ 1]*sh_b[ 2] + sh_a[ 2]*sh_b[ 4];
  c[ 5*cN+n] = sh_a[ 0]*sh_b[ 3] + sh_a[ 1]*sh_b[ 4] + sh_a[ 2]*sh_b[ 5];

  j = 3; sh_a[0] = a[n + j*aN];
  j = 4; sh_a[1] = a[n + j*aN];
  j = 5; sh_a[2] = a[n + j*aN];

  c[ 6*cN+n] = sh_a[ 0]*sh_b[ 0] + sh_a[ 1]*sh_b[ 1] + sh_a[ 2]*sh_b[ 3];
  c[ 7*cN+n] = sh_a[ 0]*sh_b[ 1] + sh_a[ 1]*sh_b[ 2] + sh_a[ 2]*sh_b[ 4];
  c[ 8*cN+n] = sh_a[ 0]*sh_b[ 3] + sh_a[ 1]*sh_b[ 4] + sh_a[ 2]*sh_b[ 5];

  j = 6; sh_a[0] = a[n + j*aN];
  j = 7; sh_a[1] = a[n + j*aN];
  j = 8; sh_a[2] = a[n + j*aN];

  c[ 9*cN+n] = sh_a[ 0]*sh_b[ 0] + sh_a[ 1]*sh_b[ 1] + sh_a[ 2]*sh_b[ 3];
  c[10*cN+n] = sh_a[ 0]*sh_b[ 1] + sh_a[ 1]*sh_b[ 2] + sh_a[ 2]*sh_b[ 4];
  c[11*cN+n] = sh_a[ 0]*sh_b[ 3] + sh_a[ 1]*sh_b[ 4] + sh_a[ 2]*sh_b[ 5];

  j = 10; sh_a[0] = a[n + j*aN];
  j = 11; sh_a[1] = a[n + j*aN];
  j = 12; sh_a[2] = a[n + j*aN];

  c[12*cN+n] = sh_a[ 0]*sh_b[ 0] + sh_a[ 1]*sh_b[ 1] + sh_a[ 2]*sh_b[ 3];
  c[13*cN+n] = sh_a[ 0]*sh_b[ 1] + sh_a[ 1]*sh_b[ 2] + sh_a[ 2]*sh_b[ 4];
  c[14*cN+n] = sh_a[ 0]*sh_b[ 3] + sh_a[ 1]*sh_b[ 4] + sh_a[ 2]*sh_b[ 5];

  j = 15; sh_a[0] = a[n + j*aN];
  j = 16; sh_a[1] = a[n + j*aN];
  j = 17; sh_a[2] = a[n + j*aN];

  c[15*cN+n] = sh_a[ 0]*sh_b[ 0] + sh_a[ 1]*sh_b[ 1] + sh_a[ 2]*sh_b[ 3];
  c[16*cN+n] = sh_a[ 0]*sh_b[ 1] + sh_a[ 1]*sh_b[ 2] + sh_a[ 2]*sh_b[ 4];
  c[17*cN+n] = sh_a[ 0]*sh_b[ 3] + sh_a[ 1]*sh_b[ 4] + sh_a[ 2]*sh_b[ 5];
}

/// Invert Cramer Symetric ////////////////////////////////////////////////////
__device__ void invertCramerSym_fn(float *a) {
  typedef float TT;

  const TT c00 = a[2] * a[5] - a[4] * a[4];
  const TT c01 = a[4] * a[3] - a[1] * a[5];
  const TT c02 = a[1] * a[4] - a[2] * a[3];
  const TT c11 = a[5] * a[0] - a[3] * a[3];
  const TT c12 = a[3] * a[1] - a[4] * a[0];
  const TT c22 = a[0] * a[2] - a[1] * a[1];

  const TT det = a[0] * c00 + a[1] * c01 + a[3] * c02;

  const TT s = TT(1) / det;

  a[0] = s*c00;
  a[1] = s*c01;
  a[2] = s*c11;
  a[3] = s*c02;
  a[4] = s*c12;
  a[5] = s*c22;
}

/// AddIntoUpperLeft3x3  //////////////////////////////////////////////////////
__device__ void addIntoUpperLeft3x3_fn(const float* __restrict__ a, size_t aN, 
                                       const float* __restrict__ b, size_t bN, 
                                       float *c, const int N) {
  int n = threadIdx.x + blockIdx.x * blockDim.x;
  
  if(n < N) {
    c[0] = a[0*aN+n] + b[0*bN+n];
    c[1] = a[1*aN+n] + b[1*bN+n];
    c[2] = a[2*aN+n] + b[2*bN+n];
    c[3] = a[3*aN+n] + b[3*bN+n];
    c[4] = a[4*aN+n] + b[4*bN+n];
    c[5] = a[5*aN+n] + b[5*bN+n];
  }
}


/// MultResidualsAdd //////////////////////////////////////////////////////////
__device__ void multResidualsAdd_fn(
    const float* __restrict__ a, size_t aN, 
    const float* __restrict__ b, size_t bN,
    const float* __restrict__ c, size_t cN,
    float *d, size_t dN, int N) {
  /*int i = threadIdx.x;*/
  int n = threadIdx.x + blockIdx.x * blockDim.x;

  float x0, x1, x2;

  /*for (int z = 0; z < (N-1)/gridDim.x  +1; z++) {*/
    /*n += z*gridDim.x;*/
     if (n < N) {
       // manually subrtact into local vars -- 3 of them
       x0 = c[0 * cN + n] - b[0 * bN + n];
       x1 = c[1 * cN + n] - b[1 * bN + n];
       x2 = c[2 * cN + n] - b[2 * bN + n];

       // generate loop (can also write it manually this time, it's not much)
       d[0 * dN + n] = b[0 * bN + n] + a[ 0 * aN + n] * x0 + a[ 1 * aN + n] * x1 + a[ 2 * aN + n] * x2;
       d[1 * dN + n] = b[1 * bN + n] + a[ 3 * aN + n] * x0 + a[ 4 * aN + n] * x1 + a[ 5 * aN + n] * x2;
       d[2 * dN + n] = b[2 * bN + n] + a[ 6 * aN + n] * x0 + a[ 7 * aN + n] * x1 + a[ 8 * aN + n] * x2;
       d[3 * dN + n] = b[3 * bN + n] + a[ 9 * aN + n] * x0 + a[10 * aN + n] * x1 + a[11 * aN + n] * x2;
       d[4 * dN + n] = b[4 * bN + n] + a[12 * aN + n] * x0 + a[13 * aN + n] * x1 + a[14 * aN + n] * x2;
       d[5 * dN + n] = b[5 * bN + n] + a[15 * aN + n] * x0 + a[16 * aN + n] * x1 + a[17 * aN + n] * x2;
     }
  /*}*/
}

/// KalmanGain_x_propErr //////////////////////////////////////////////////////
__device__ void kalmanGain_x_propErr_fn(
    const float* __restrict__ d_kalmanGain, size_t stride_kalmanGain,
    const float* __restrict__ d_propErr, size_t stride_propErr,
    float *d_outErr, size_t stride_outErr, const int N) {
  // a = d_kalmanGain,  b = d_propErr, c = outErrTemp
  // c = b - a*b
  const float *a = d_kalmanGain;
  const float *b = d_propErr;
  float *c = d_outErr;

  size_t aN = stride_kalmanGain;
  size_t bN = stride_propErr;
  size_t cN = stride_outErr;

  __shared__ float sh_a[LH][BLOCK_SIZE_X];
  __shared__ float sh_b[LS][BLOCK_SIZE_X];
  
  register float reg_c[LS];

  int i = threadIdx.x;
  int n = threadIdx.x + blockIdx.x * blockDim.x;

  if (n < N) {
      for (int j = 0; j < LH; j++) {
        sh_a[j][i] = a[n + j*aN];
      }
      for (int j = 0; j < LS; j++) {
        sh_b[j][i] = b[n + j*bN];
      }
      __syncthreads();
      reg_c[ 0] = sh_a[ 0][i]*sh_b[ 0][i] + sh_a[ 1][i]*sh_b[ 1][i] + sh_a[ 2][i]*sh_b[ 3][i];
      reg_c[ 1] = sh_a[ 3][i]*sh_b[ 0][i] + sh_a[ 4][i]*sh_b[ 1][i] + sh_a[ 5][i]*sh_b[ 3][i];
      reg_c[ 2] = sh_a[ 3][i]*sh_b[ 1][i] + sh_a[ 4][i]*sh_b[ 2][i] + sh_a[ 5][i]*sh_b[ 4][i];
      reg_c[ 3] = sh_a[ 6][i]*sh_b[ 0][i] + sh_a[ 7][i]*sh_b[ 1][i] + sh_a[ 8][i]*sh_b[ 3][i];
      reg_c[ 4] = sh_a[ 6][i]*sh_b[ 1][i] + sh_a[ 7][i]*sh_b[ 2][i] + sh_a[ 8][i]*sh_b[ 4][i];
      reg_c[ 5] = sh_a[ 6][i]*sh_b[ 3][i] + sh_a[ 7][i]*sh_b[ 4][i] + sh_a[ 8][i]*sh_b[ 5][i];
      reg_c[ 6] = sh_a[ 9][i]*sh_b[ 0][i] + sh_a[10][i]*sh_b[ 1][i] + sh_a[11][i]*sh_b[ 3][i];
      reg_c[ 7] = sh_a[ 9][i]*sh_b[ 1][i] + sh_a[10][i]*sh_b[ 2][i] + sh_a[11][i]*sh_b[ 4][i];
      reg_c[ 8] = sh_a[ 9][i]*sh_b[ 3][i] + sh_a[10][i]*sh_b[ 4][i] + sh_a[11][i]*sh_b[ 5][i];
      reg_c[ 9] = sh_a[ 9][i]*sh_b[ 6][i] + sh_a[10][i]*sh_b[ 7][i] + sh_a[11][i]*sh_b[ 8][i];
      reg_c[10] = sh_a[12][i]*sh_b[ 0][i] + sh_a[13][i]*sh_b[ 1][i] + sh_a[14][i]*sh_b[ 3][i];
      reg_c[11] = sh_a[12][i]*sh_b[ 1][i] + sh_a[13][i]*sh_b[ 2][i] + sh_a[14][i]*sh_b[ 4][i];
      reg_c[12] = sh_a[12][i]*sh_b[ 3][i] + sh_a[13][i]*sh_b[ 4][i] + sh_a[14][i]*sh_b[ 5][i];
      reg_c[13] = sh_a[12][i]*sh_b[ 6][i] + sh_a[13][i]*sh_b[ 7][i] + sh_a[14][i]*sh_b[ 8][i];
      reg_c[14] = sh_a[12][i]*sh_b[10][i] + sh_a[13][i]*sh_b[11][i] + sh_a[14][i]*sh_b[12][i];
      reg_c[15] = sh_a[15][i]*sh_b[ 0][i] + sh_a[16][i]*sh_b[ 1][i] + sh_a[17][i]*sh_b[ 3][i];
      reg_c[16] = sh_a[15][i]*sh_b[ 1][i] + sh_a[16][i]*sh_b[ 2][i] + sh_a[17][i]*sh_b[ 4][i];
      reg_c[17] = sh_a[15][i]*sh_b[ 3][i] + sh_a[16][i]*sh_b[ 4][i] + sh_a[17][i]*sh_b[ 5][i];
      reg_c[18] = sh_a[15][i]*sh_b[ 6][i] + sh_a[16][i]*sh_b[ 7][i] + sh_a[17][i]*sh_b[ 8][i];
      reg_c[19] = sh_a[15][i]*sh_b[10][i] + sh_a[16][i]*sh_b[11][i] + sh_a[17][i]*sh_b[12][i];
      reg_c[20] = sh_a[15][i]*sh_b[15][i] + sh_a[16][i]*sh_b[16][i] + sh_a[17][i]*sh_b[17][i];
      // TODO: profile to see if reg_c values are really in registers.
#pragma unroll
      for (int j = 0; j < LS; j++) {
        /*d_outErr[j*N + n] = d_propErr[j*N + n] - d_outErr[j*N + n];*/
        c[j*cN + n] = sh_b[j][i] - reg_c[j];
      }
   }
}

__global__ void kalmanUpdateMerged_kernel(const float* __restrict__ propErr, size_t propErr_stride,
                                          const float* __restrict__ msErr, size_t msErr_stride,
                                          float *kalmanGain, size_t kalmanGain_stride,
                                          const float* __restrict__ par_iP, size_t par_iP_stride,
                                          const float* __restrict__ msPar, size_t msPar_stride,
                                          float *par_iC, size_t par_iC_stride,
                                          float *outErr, size_t outErr_stride,
                                          const int N) {
  int n = threadIdx.x + blockIdx.x * blockDim.x;
  float resErr_reg[HS];

  if (n < N) {
    for (int j = 0; j < HS; ++j) {
      resErr_reg[j] = 0; //resErr[j*resErr_stride + n];
    }
    addIntoUpperLeft3x3_fn(propErr, propErr_stride,
        msErr, msErr_stride, resErr_reg, N);
    invertCramerSym_fn(resErr_reg);

    upParam_MultKalmanGain_raw(
        propErr, propErr_stride, resErr_reg, 0,
        kalmanGain, kalmanGain_stride, N);             

    multResidualsAdd_fn(kalmanGain, kalmanGain_stride, par_iP, par_iP_stride, 
        msPar, msPar_stride, par_iC, par_iC_stride, N);

    kalmanGain_x_propErr_fn(kalmanGain, kalmanGain_stride,
                            propErr, propErr_stride,
                            outErr, outErr_stride, N);
  }
}

void kalmanUpdateMerged_wrapper(cudaStream_t& stream,
    GPlex<float>& d_propErr, GPlex<float>& d_msErr,
    GPlex<float>& d_kalmanGain, 
    GPlex<float>& d_par_iP, GPlex<float>& d_msPar,
    GPlex<float>& d_par_iC, GPlex<float>& d_outErr,
    const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  kalmanUpdateMerged_kernel <<<grid, block, 0, stream >>>
      (d_propErr.ptr, d_propErr.stride,
       d_msErr.ptr, d_msErr.stride,
       d_kalmanGain.ptr, d_kalmanGain.stride,
       d_par_iP.ptr, d_par_iP.stride,
       d_msPar.ptr, d_msPar.stride, 
       d_par_iC.ptr, d_par_iC.stride, 
       d_outErr.ptr, d_outErr.stride,
       N);
}
