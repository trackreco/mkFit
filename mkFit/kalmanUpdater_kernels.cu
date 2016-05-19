#include "Config.h"
#include "kalmanUpdater_kernels.h"

// TODO: Clean all the hard-coded #define
#define LS 21
#define HS 6
#define LH 18

#define BLOCK_SIZE_X 32
#define MAX_BLOCKS_X 65535 // CUDA constraint


/// MultKalmanGain ////////////////////////////////////////////////////////////

__device__ void upParam_MultKalmanGain_fn(
    const float* __restrict__ a, size_t aN,
    float* b_reg, float *c, int N, int n) {
  // const T* __restrict__ tells the compiler that it can uses the read-only
  // cache, without worrying about coherency.
  // c -> kalmanGain, in register

  /*int n = threadIdx.x + blockIdx.x * blockDim.x;*/
  // use registers to store values of 'a' that are use multiple times
  // To reduce the number of registers and avoid spilling interlace
  // read and compute instructions.
  float a_reg[3];

  int j;
  j = 0; a_reg[0] = a[n + j*aN];
  j = 1; a_reg[1] = a[n + j*aN];
  j = 3; a_reg[2] = a[n + j*aN];

  c[ 0] = a_reg[ 0]*b_reg[ 0] + a_reg[ 1]*b_reg[ 1] + a_reg[ 2]*b_reg[ 3];
  c[ 1] = a_reg[ 0]*b_reg[ 1] + a_reg[ 1]*b_reg[ 2] + a_reg[ 2]*b_reg[ 4];
  c[ 2] = a_reg[ 0]*b_reg[ 3] + a_reg[ 1]*b_reg[ 4] + a_reg[ 2]*b_reg[ 5];

  j = 1; a_reg[0] = a[n + j*aN];
  j = 2; a_reg[1] = a[n + j*aN];
  j = 4; a_reg[2] = a[n + j*aN];

  c[ 3] = a_reg[ 0]*b_reg[ 0] + a_reg[ 1]*b_reg[ 1] + a_reg[ 2]*b_reg[ 3];
  c[ 4] = a_reg[ 0]*b_reg[ 1] + a_reg[ 1]*b_reg[ 2] + a_reg[ 2]*b_reg[ 4];
  c[ 5] = a_reg[ 0]*b_reg[ 3] + a_reg[ 1]*b_reg[ 4] + a_reg[ 2]*b_reg[ 5];

  j = 3; a_reg[0] = a[n + j*aN];
  j = 4; a_reg[1] = a[n + j*aN];
  j = 5; a_reg[2] = a[n + j*aN];

  c[ 6] = a_reg[ 0]*b_reg[ 0] + a_reg[ 1]*b_reg[ 1] + a_reg[ 2]*b_reg[ 3];
  c[ 7] = a_reg[ 0]*b_reg[ 1] + a_reg[ 1]*b_reg[ 2] + a_reg[ 2]*b_reg[ 4];
  c[ 8] = a_reg[ 0]*b_reg[ 3] + a_reg[ 1]*b_reg[ 4] + a_reg[ 2]*b_reg[ 5];

  j = 6; a_reg[0] = a[n + j*aN];
  j = 7; a_reg[1] = a[n + j*aN];
  j = 8; a_reg[2] = a[n + j*aN];

  c[ 9] = a_reg[ 0]*b_reg[ 0] + a_reg[ 1]*b_reg[ 1] + a_reg[ 2]*b_reg[ 3];
  c[10] = a_reg[ 0]*b_reg[ 1] + a_reg[ 1]*b_reg[ 2] + a_reg[ 2]*b_reg[ 4];
  c[11] = a_reg[ 0]*b_reg[ 3] + a_reg[ 1]*b_reg[ 4] + a_reg[ 2]*b_reg[ 5];

  j = 10; a_reg[0] = a[n + j*aN];
  j = 11; a_reg[1] = a[n + j*aN];
  j = 12; a_reg[2] = a[n + j*aN];

  c[12] = a_reg[ 0]*b_reg[ 0] + a_reg[ 1]*b_reg[ 1] + a_reg[ 2]*b_reg[ 3];
  c[13] = a_reg[ 0]*b_reg[ 1] + a_reg[ 1]*b_reg[ 2] + a_reg[ 2]*b_reg[ 4];
  c[14] = a_reg[ 0]*b_reg[ 3] + a_reg[ 1]*b_reg[ 4] + a_reg[ 2]*b_reg[ 5];

  j = 15; a_reg[0] = a[n + j*aN];
  j = 16; a_reg[1] = a[n + j*aN];
  j = 17; a_reg[2] = a[n + j*aN];

  c[15] = a_reg[ 0]*b_reg[ 0] + a_reg[ 1]*b_reg[ 1] + a_reg[ 2]*b_reg[ 3];
  c[16] = a_reg[ 0]*b_reg[ 1] + a_reg[ 1]*b_reg[ 2] + a_reg[ 2]*b_reg[ 4];
  c[17] = a_reg[ 0]*b_reg[ 3] + a_reg[ 1]*b_reg[ 4] + a_reg[ 2]*b_reg[ 5];
}

/// Invert Cramer Symetric ////////////////////////////////////////////////////
__device__ void invertCramerSym_fn(float *a) {
  // a is in registers.
  // to use global memory, a stride will be required and accesses would be:
  // a[n + stride_a * i];
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
                                       float *c, const int N, int n) {
  /*int n = threadIdx.x + blockIdx.x * blockDim.x;*/
  
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
    float* reg_a,
    const float* __restrict__ b, size_t bN,
    const float* __restrict__ c, size_t cN,
    float *d, size_t dN, int N, int n) {
  // a -> kalmanGain

  /*int i = threadIdx.x;*/
  /*int n = threadIdx.x + blockIdx.x * blockDim.x;*/

  /*for (int z = 0; z < (N-1)/gridDim.x  +1; z++) {*/
    /*n += z*gridDim.x;*/
     if (n < N) {
       // manually substract into local vars -- 3 of them
       const float x0 = c[0 * cN + n] - b[0 * bN + n];
       const float x1 = c[1 * cN + n] - b[1 * bN + n];
       const float x2 = c[2 * cN + n] - b[2 * bN + n];

       // generate loop (can also write it manually this time, it's not much)
       // WARNING: highly numerically sensitive expressions.
       //          order of operand yield different results.
       //          probably also true for x0, x1, x2.
       //          check if that has an impact on the physics
       //          
       //          use of register allow ~35% better performance; but
       //          it does have a small impact on pT's mean and std deviation
       d[0 * dN + n] = b[0 * bN + n] + (reg_a[ 0] * x0 + reg_a[ 1] * x1 + reg_a[ 2] * x2);
       d[1 * dN + n] = b[1 * bN + n] + (reg_a[ 3] * x0 + reg_a[ 4] * x1 + reg_a[ 5] * x2);
       d[2 * dN + n] = b[2 * bN + n] + (reg_a[ 6] * x0 + reg_a[ 7] * x1 + reg_a[ 8] * x2);
       d[3 * dN + n] = b[3 * bN + n] + (reg_a[ 9] * x0 + reg_a[10] * x1 + reg_a[11] * x2);
       d[4 * dN + n] = b[4 * bN + n] + (reg_a[12] * x0 + reg_a[13] * x1 + reg_a[14] * x2);
       d[5 * dN + n] = b[5 * bN + n] + (reg_a[15] * x0 + reg_a[16] * x1 + reg_a[17] * x2);
     }
  /*}*/
}

/// KalmanGain_x_propErr //////////////////////////////////////////////////////
__device__ void kalmanGain_x_propErr_fn(
    float* d_kalmanGain,
    const float* __restrict__ d_propErr, size_t stride_propErr,
    float *d_outErr, size_t stride_outErr, const int N, int n) {
  // a = d_kalmanGain,  b = d_propErr, c = outErrTemp
  // c = b - a*b
  float *a = d_kalmanGain;
  const float *b = d_propErr;
  float *c = d_outErr;

  /*size_t aN = stride_kalmanGain;*/
  size_t bN = stride_propErr;
  size_t cN = stride_outErr;

  register float reg_c[LS];

  /*int i = threadIdx.x;*/
  /*int n = threadIdx.x + blockIdx.x * blockDim.x;*/

  if (n < N) {
      // Warning: parenthesis are required to ensure that
      //          the result is the same than with two separate steps.
      // Note: It seems slightly faster to explicitly keep c in registers
      //       and fill c with them later on.
      reg_c[ 0] =  b[n+  0*bN] - (a[ 0]*b[n +  0*bN] + a[ 1]*b[n +  1*bN] + a[ 2]*b[n +  3*bN]);
      reg_c[ 1] =  b[n+  1*bN] - (a[ 3]*b[n +  0*bN] + a[ 4]*b[n +  1*bN] + a[ 5]*b[n +  3*bN]);
      reg_c[ 2] =  b[n+  2*bN] - (a[ 3]*b[n +  1*bN] + a[ 4]*b[n +  2*bN] + a[ 5]*b[n +  4*bN]);
      reg_c[ 3] =  b[n+  3*bN] - (a[ 6]*b[n +  0*bN] + a[ 7]*b[n +  1*bN] + a[ 8]*b[n +  3*bN]);
      reg_c[ 4] =  b[n+  4*bN] - (a[ 6]*b[n +  1*bN] + a[ 7]*b[n +  2*bN] + a[ 8]*b[n +  4*bN]);
      reg_c[ 5] =  b[n+  5*bN] - (a[ 6]*b[n +  3*bN] + a[ 7]*b[n +  4*bN] + a[ 8]*b[n +  5*bN]);
      reg_c[ 6] =  b[n+  6*bN] - (a[ 9]*b[n +  0*bN] + a[10]*b[n +  1*bN] + a[11]*b[n +  3*bN]);
      reg_c[ 7] =  b[n+  7*bN] - (a[ 9]*b[n +  1*bN] + a[10]*b[n +  2*bN] + a[11]*b[n +  4*bN]);
      reg_c[ 8] =  b[n+  8*bN] - (a[ 9]*b[n +  3*bN] + a[10]*b[n +  4*bN] + a[11]*b[n +  5*bN]);
      reg_c[ 9] =  b[n+  9*bN] - (a[ 9]*b[n +  6*bN] + a[10]*b[n +  7*bN] + a[11]*b[n +  8*bN]);
      reg_c[10] =  b[n+ 10*bN] - (a[12]*b[n +  0*bN] + a[13]*b[n +  1*bN] + a[14]*b[n +  3*bN]);
      reg_c[11] =  b[n+ 11*bN] - (a[12]*b[n +  1*bN] + a[13]*b[n +  2*bN] + a[14]*b[n +  4*bN]);
      reg_c[12] =  b[n+ 12*bN] - (a[12]*b[n +  3*bN] + a[13]*b[n +  4*bN] + a[14]*b[n +  5*bN]);
      reg_c[13] =  b[n+ 13*bN] - (a[12]*b[n +  6*bN] + a[13]*b[n +  7*bN] + a[14]*b[n +  8*bN]);
      reg_c[14] =  b[n+ 14*bN] - (a[12]*b[n + 10*bN] + a[13]*b[n + 11*bN] + a[14]*b[n + 12*bN]);
      reg_c[15] =  b[n+ 15*bN] - (a[15]*b[n +  0*bN] + a[16]*b[n +  1*bN] + a[17]*b[n +  3*bN]);
      reg_c[16] =  b[n+ 16*bN] - (a[15]*b[n +  1*bN] + a[16]*b[n +  2*bN] + a[17]*b[n +  4*bN]);
      reg_c[17] =  b[n+ 17*bN] - (a[15]*b[n +  3*bN] + a[16]*b[n +  4*bN] + a[17]*b[n +  5*bN]);
      reg_c[18] =  b[n+ 18*bN] - (a[15]*b[n +  6*bN] + a[16]*b[n +  7*bN] + a[17]*b[n +  8*bN]);
      reg_c[19] =  b[n+ 19*bN] - (a[15]*b[n + 10*bN] + a[16]*b[n + 11*bN] + a[17]*b[n + 12*bN]);
      reg_c[20] =  b[n+ 20*bN] - (a[15]*b[n + 15*bN] + a[16]*b[n + 16*bN] + a[17]*b[n + 17*bN]);
#pragma unroll
      for (int j = 0; j < LS; j++) {
        c[j*cN + n] = reg_c[j];
      }
   }
}

__global__ void kalmanUpdate_kernel(
    const float* __restrict__ propErr, size_t propErr_stride,
    const float* __restrict__ msErr, size_t msErr_stride,
    const float* __restrict__ par_iP, size_t par_iP_stride,
    const float* __restrict__ msPar, size_t msPar_stride,
    float *par_iC, size_t par_iC_stride,
    float *outErr, size_t outErr_stride,
    const int N) {
  int grid_width = blockDim.x * gridDim.x;
  // Note: similar results with propErr kept in registers.
  //       It is read-only so using the read-only cache yields more flexibility
  //       wrt block size without increasing the pressure on registers to much.
  int n = threadIdx.x + blockIdx.x * blockDim.x;
  // There is no need to keep resErr and kalmanGain as global memory arrays.
  float resErr_reg[HS];
  float kalmanGain_reg[LH];

  // If there is more matrices than MAX_BLOCKS_X * BLOCK_SIZE_X 
  for (int z = 0; z < (N-1)/grid_width  +1; z++) {
    /*n += z*gridDim.x;*/
    n += z*grid_width;
    if (n < N) {
      for (int j = 0; j < HS; ++j) {
        resErr_reg[j] = 0; //resErr[j*resErr_stride + n];
      }
      addIntoUpperLeft3x3_fn(propErr, propErr_stride,
          msErr, msErr_stride, resErr_reg, N, n);
      invertCramerSym_fn(resErr_reg);

      upParam_MultKalmanGain_fn(propErr, propErr_stride,
          resErr_reg, kalmanGain_reg, N, n);             
      multResidualsAdd_fn(kalmanGain_reg, par_iP, par_iP_stride, 
          msPar, msPar_stride, par_iC, par_iC_stride, N, n);

      kalmanGain_x_propErr_fn(kalmanGain_reg,
          propErr, propErr_stride,
          outErr, outErr_stride, N, n);
    }
  }
}

void kalmanUpdate_wrapper(cudaStream_t& stream,
    GPlex<float>& d_propErr, GPlex<float>& d_msErr,
    GPlex<float>& d_par_iP, GPlex<float>& d_msPar,
    GPlex<float>& d_par_iC, GPlex<float>& d_outErr,
    const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  kalmanUpdate_kernel <<<grid, block, 0, stream >>>
      (d_propErr.ptr, d_propErr.stride,
       d_msErr.ptr, d_msErr.stride,
       d_par_iP.ptr, d_par_iP.stride,
       d_msPar.ptr, d_msPar.stride, 
       d_par_iC.ptr, d_par_iC.stride, 
       d_outErr.ptr, d_outErr.stride,
       N);
}

// Should probably not be in this file, but creating a file for
// this oneliner seems overkill.
void separate_first_call_for_meaningful_profiling_numbers() {
  cudaDeviceSynchronize();
}
