#include "Config.h"
#include "Hit.h"
#include "kalmanUpdater_kernels.h"
#include "computeChi2_kernels.h"
#include "gpu_utils.h"

// TODO: Clean all the hard-coded #define
#define LS 21
#define HS 6
#define LH 18
#define HV 3

#define BLOCK_SIZE_X 32


__device__ void subtract_matrix(const float *a, const int aN, 
                                const float *b, const int bN, 
                                      float *c, const int cN,
                                const int size, const int n) {
  for (int i = 0; i < size; ++i) {
    c[i*cN + n] = a[i*aN + n] - b[i*bN + n];
    
  }
}

__device__ float getHypot_fn(const float x, const float y)
{
  return sqrt(x*x + y*y);
}

__device__
void KalmanHTG_fn(const GPlexRegQF& a00, const GPlexRegQF& a01,
                  const GPlexReg2S &b, GPlexRegHH &c)
{

   // HTG  = rot * res_loc
   //   C  =  A  *    B   

   // Based on script generation and adapted to custom sizes.
      c[ 0] = a00[0]*b[ 0];
      c[ 1] = a00[0]*b[ 1];
      c[ 2] = 0.;
      c[ 3] = a01[0]*b[ 0];
      c[ 4] = a01[0]*b[ 1];
      c[ 5] = 0.;
      c[ 6] = b[ 1];
      c[ 7] = b[ 2];
      c[ 8] = 0.;
}

__device__
void KalmanGain_fn(const GPlexLS &A, const GPlexRegHH &b, GPlexRegLH &c, const int n)
{
  // C = A * B, C is 6x3, A is 6x6 sym , B is 6x3
  using T = float;
  float *a = A.ptr;
  int aN = A.stride; int an = n;  // Global array
  int bN = 1;        int bn = 0;  // Register array
  int cN = 1;        int cn = 0;

#include "KalmanGain.ah"
}

#include "KalmanUtilsMPlex.icc"

__device__
void KHMult_fn(const GPlexRegLH &A, 
               const GPlexRegQF& B00,
               const GPlexRegQF& B01,
               GPlexRegLL &C)
{
  KHMult_imp(A, B00, B01, C, 0, 1);
}

__device__
void KHC_fn(const GPlexRegLL &a, const GPlexLS &B, GPlexLS &C, const int n)
{
  // C = A * B, C is 6x6, A is 6x6 , B is 6x6 sym
  using T = float;
                 int aN = 1; int an = 0;  // Register array
  T *b = B.ptr;  int bN = B.stride;  int bn = n;
  T *c = C.ptr;  int cN = C.stride;  int cn = n;
#include "KHC.ah"
}

// 
__device__
void ConvertToCCS_fn(const GPlexLV &a, GPlexRegLV &b, GPlexRegLL &c, const int n)
{
  ConvertToCCS_imp(a, b, c, n, n+1);
}

__device__
void PolarErr_fn(const GPlexRegLL &a, const float *b, int bN, GPlexRegLL &c, const int n)
{
  // C = A * B, C is 6x6, A is 6x6 , B is 6x6 sym
 
  // Generated code access arrays with variables cN, cn
  // c[i*cN+cn]  
  int aN = 1; int an = 0;  // Register array
              int bn = n;  // Global array
  int cN = 1; int cn = 0;
#include "CCSErr.ah"
}

__device__
void PolarErrTransp_fn(const GPlexRegLL &a, const GPlexRegLL &b, GPlexLS &C, const int n)
{
  // C = A * B, C is sym, A is 6x6 , B is 6x6
  using T = float;
                 int aN = 1;         int an = 0;
                 int bN = 1;         int bn = 0;
  T *c = C.ptr;  int cN = C.stride;  int cn = n;
#include "CCSErrTransp.ah"
}

__device__
void ConvertToCartesian_fn(const GPlexRegLV &a, GPlexLV& b, GPlexRegLL &c, const int n)
{
  ConvertToCartesian_imp(a, b, c, n, n+1); 
}

__device__
void CartesianErr_fn(const GPlexRegLL &a, const float *b, const int bN, GPlexRegLL &c, const int n)
{
  // C = A * B, C is 6x6, A is 6x6 , B is 6x6 sym
  int aN = 1; int an = 0;
              int bn = n;
  int cN = 1; int cn = 0;

#include "CartesianErr.ah"
}

__device__
void CartesianErrTransp_fn(const GPlexRegLL &a, const GPlexRegLL &b, GPlexLS &C, const int n)
{
  // C = A * B, C is sym, A is 6x6 , B is 6x6
  using T = float;
  int aN = 1; int an = 0;
  int bN = 1; int bn = 0;
  T *c = C.ptr;  int cN = C.stride;  int cn = n;

#include "CartesianErrTransp.ah"
}


/// MultKalmanGain ////////////////////////////////////////////////////////////

__device__ void upParam_MultKalmanGain_fn(
    const float* __restrict__ a, const size_t aN,
    const float* b_reg, float *c, const int N, const int n) {
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

__device__ void invertCramerSym2x2_fn(GPlexReg2S &a) {
  float det = a[0] * a[2] - a[1] * a[1];
  const float s   = float(1) / det;
  const float tmp = s * a[2];
  a[1] *= -s;
  a[2]  = s * a[0];
  a[0]  = tmp;
}

__device__ void subtractFirst3_fn(const GPlexHV __restrict__ &A,
                                  const GPlexLV __restrict__ &B,
                                  GPlexRegHV &C, const int N, int n) {
  using T = float;
  const T *a = A.ptr;  int aN = A.stride;
  const T *b = B.ptr;  int bN = B.stride;
        T *c = C.arr;
  /*int n = threadIdx.x + blockIdx.x * blockDim.x;*/
  
  if(n < N) {
    c[0] = a[0*aN+n] - b[0*bN+n];
    c[1] = a[1*aN+n] - b[1*bN+n];
    c[2] = a[2*aN+n] - b[2*bN+n];
  }
}

/// AddIntoUpperLeft3x3  //////////////////////////////////////////////////////
__device__ void addIntoUpperLeft3x3_fn(const GPlexLS __restrict__ &A,
                                       const GPlexHS __restrict__ &B,
                                       GPlexRegHS &C, const int N, const int n) {
  using T = float;
  const T *a = A.ptr;  int aN = A.stride;
  const T *b = B.ptr;  int bN = B.stride;
        T *c = C.arr;
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
    const GPlexRegLH &reg_a,
    const GPlexLV __restrict__ &B, 
    const GPlexReg2V &c,
          GPlexLV &D,
    const int N, const int n) {

  MultResidualsAdd_imp(reg_a, B, c, D, n, n+1);
}

__device__
void MultResidualsAdd_all_reg(const GPlexRegLH &a,
		      const GPlexRegLV &b,
		      const GPlexReg2V &c,
          GPlexRegLV &d)
{
   // outPar = psPar + kalmanGain*(dPar)
   //   D    =   B         A         C
   // where right half of kalman gain is 0 

   // XXX Regenerate with a script.
      // generate loop (can also write it manually this time, it's not much)
      d[0] = b[0] + a[ 0] * c[0] + a[ 1] * c[1];
      d[1] = b[1] + a[ 3] * c[0] + a[ 4] * c[1];
      d[2] = b[2] + a[ 6] * c[0] + a[ 7] * c[1];
      d[3] = b[3] + a[ 9] * c[0] + a[10] * c[1];
      d[4] = b[4] + a[12] * c[0] + a[13] * c[1];
      d[5] = b[5] + a[15] * c[0] + a[16] * c[1];
}


__device__ void kalmanUpdate_fn(
    GPlexLS &propErr, const GPlexHS __restrict__ &msErr,
    const GPlexLV __restrict__ &par_iP, const GPlexHV __restrict__ &msPar,
    GPlexLV &par_iC, GPlexLS &outErr, const int n, const int N) {
  // Note: similar results with propErr kept in registers.
  //       It is read-only so using the read-only cache yields more flexibility
  //       wrt block size without increasing the pressure on registers to much.
  // There is no need to keep resErr and kalmanGain as global memory arrays.
  GPlexRegHS resErr_reg;

  // If there is more matrices than max_blocks_x * BLOCK_SIZE_X 
  if (n < N) {
    for (int j = 0; j < HS; ++j) {
      resErr_reg[j] = 0; //resErr[j*resErr_stride + n];
    }

    // FIXME: Add useCMSGeom -> port propagateHelixToRMPlex
#if 0
    if (Config::useCMSGeom) {
      propagateHelixToRMPlex(psErr,  psPar, inChg,  msPar, propErr, propPar);
    } else {
      propErr = psErr;
      propPar = psPar;
    }
#endif
    GPlexRegQF rotT00;
    GPlexRegQF rotT01;
    const float r = hipo(msPar(n, 0, 0), msPar(n, 1, 0));
    rotT00[0] = -(msPar(n, 1, 0) + par_iP(n, 1, 0))/(2*r);
    rotT01[0] =  (msPar(n, 0, 0) + par_iP(n, 0, 0))/(2*r);

    GPlexRegHV res_glo;
    subtractFirst3_fn(msPar, par_iP, res_glo, N, n);

    addIntoUpperLeft3x3_fn(propErr, msErr, resErr_reg, N, n);
    GPlexReg2V res_loc;   //position residual in local coordinates
    RotateResidulsOnTangentPlane_fn(rotT00,rotT01,res_glo,res_loc);
    GPlexReg2S resErr_loc; // 2x2 sym
    GPlexRegHH tempHH;  // 3*3 sym
    ProjectResErr_fn  (rotT00, rotT01, resErr_reg, tempHH);
    ProjectResErrTransp_fn(rotT00, rotT01, tempHH, resErr_loc);

    /*invertCramerSym_fn(resErr_reg);*/
    invertCramerSym2x2_fn(resErr_loc);
#ifndef CCSCOORD
    // Move to "polar" coordinates: (x,y,z,1/pT,phi,theta) [can we find a better name?]

    GPlexRegLV propPar_pol;  // propagated parameters in "polar" coordinates*/
    GPlexRegLL jac_pol;  // jacobian from cartesian to "polar"*/
    ConvertToCCS_fn(par_iP, propPar_pol, jac_pol, n);

    GPlexRegLL tempLL;
    PolarErr_fn(jac_pol, propErr.ptr, propErr.stride, tempLL, n);
    PolarErrTransp_fn(jac_pol, tempLL, propErr, n);// propErr is now propagated errors in "polar" coordinates
#endif

    // Kalman update in "polar" coordinates
    GPlexRegLH K;
    KalmanHTG_fn(rotT00, rotT01, resErr_loc, tempHH);
    KalmanGain_fn(propErr, tempHH, K, n);

#ifdef CCSCOORD
    multResidualsAdd_fn(K, par_iP, res_loc, par_iC, N, n);// propPar_pol is now the updated parameters in "polar" coordinates
    GPlexRegLL tempLL;
#else
    MultResidualsAdd_all_reg(K, propPar_pol, res_loc, propPar_pol);
#endif

    KHMult_fn(K, rotT00, rotT01, tempLL);
    KHC_fn(tempLL, propErr, outErr, n);
    subtract_matrix(propErr.ptr, propErr.stride, outErr.ptr, outErr.stride, 
                    outErr.ptr, outErr.stride, LS, n);

#ifndef CCSCOORD
    // Go back to cartesian coordinates

    // jac_pol is now the jacobian from "polar" to cartesian
    // outPar -> par_iC
    ConvertToCartesian_fn(propPar_pol, par_iC, jac_pol, n);
    CartesianErr_fn      (jac_pol, outErr.ptr, outErr.stride, tempLL, n);
    CartesianErrTransp_fn(jac_pol, tempLL, outErr, n);// outErr is in cartesian coordinates now
#endif
  }
}

__global__ void kalmanUpdate_kernel(
    GPlexLS propErr, const GPlexHS __restrict__ msErr,
    const GPlexLV __restrict__ par_iP, const GPlexHV __restrict__ msPar,
    GPlexLV par_iC, GPlexLS outErr, const int N) {
  int grid_width = blockDim.x * gridDim.x;
  int n = threadIdx.x + blockIdx.x * blockDim.x;

  for (int z = 0; z < (N-1)/grid_width  +1; z++) {
    n += z*grid_width;
    kalmanUpdate_fn(propErr, msErr, par_iP, msPar, par_iC, outErr, n, N);
  }
}

void kalmanUpdate_wrapper(const cudaStream_t& stream,
    GPlexLS& d_propErr, const GPlexHS& d_msErr,
    GPlexLV& d_par_iP, const GPlexHV& d_msPar,
    GPlexLV& d_par_iC, GPlexLS& d_outErr,
    const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  kalmanUpdate_kernel <<<grid, block, 0, stream >>>
      (d_propErr, d_msErr, d_par_iP, d_msPar, d_par_iC, d_outErr, N);
}

