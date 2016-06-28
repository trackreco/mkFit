#include "computeChi2_kernels.h"

#include "GPlex.h"
#include "kalmanUpdater_kernels.h"
#include "gpu_utils.h"

#include <stdio.h>

#define L 6
#define HS 6
#define HV 3
#define BLOCK_SIZE_X 256


__device__ void chi2Similarity_fn(
    const GPlexReg2V &a,
    const GPlexReg2S &c, // in registers
    float *d, const size_t dN) {

  int n = threadIdx.x + blockIdx.x * blockDim.x;

  // manually subrtact into local vars -- 3 of them
  /*float x0 = a[0 * aN + n] - b[0 * aN + n];*/
  /*float x1 = a[1 * aN + n] - b[1 * aN + n];*/
  /*float x2 = a[2 * aN + n] - b[2 * aN + n];*/
  /*d[0 * dN + n] = c[0]*x0*x0 + c[2]*x1*x1 + c[5]*x2*x2 +*/
              /*2*( c[1]*x1*x0 + c[3]*x2*x0 + c[4]*x1*x2);*/
  d[0 * dN + n] = c[0]*a[0]*a[0]
                + c[2]*a[1]*a[1] 
            + 2*( c[1]*a[1]*a[0]);
}


__device__ void RotateResidulsOnTangentPlane_fn(const float r00,//r00
				  const float r01,//r01
				  const GPlexRegHV &a  ,//res_glo
          GPlexReg2V &b  )//res_loc
{

   // res_loc = rotT * res_glo
   //   B     =  R   *    A   
  b[0] =  r00*a[0] + r01*a[1];
  b[1] =  a[2];
}


__device__ void ProjectResErr_fn(const float a00,
		   const float a01,
		   const GPlexRegHS &b, 
       GPlexRegHH &c)
{
  // C = A * B, C is 3x3, A is 3x3 , B is 3x3 sym

  // Based on script generation and adapted to custom sizes.
      c[ 0] = a00*b[ 0] + a01*b[ 1];
      c[ 1] = a00*b[ 1] + a01*b[ 2];
      c[ 2] = a00*b[ 3] + a01*b[ 4];
      c[ 3] = b[ 3];
      c[ 4] = b[ 4];
      c[ 5] = b[ 5];
      c[ 6] = a01*b[ 0] - a00*b[ 1];
      c[ 7] = a01*b[ 1] - a00*b[ 2];
      c[ 8] = a01*b[ 3] - a00*b[ 4];
}


__device__ void ProjectResErrTransp_fn(const float a00,
			 const float a01, const GPlexRegHH &b, GPlexReg2S &c)
{
  // C = A * B, C is 3x3 sym, A is 3x3 , B is 3x3

  // Based on script generation and adapted to custom sizes.
      c[ 0] = b[ 0]*a00 + b[ 1]*a01;
      c[ 1] = b[ 3]*a00 + b[ 4]*a01;
      c[ 2] = b[ 5];
}


__device__ void computeChi2_fn(
    const GPlexLS &propErr, const GPlexHS &msErr, const GPlexHV &msPar,
    const GPlexLV &propPar, GPlexQF &outChi2, const int N) {
  int grid_width = blockDim.x * gridDim.x;
  int n = threadIdx.x + blockIdx.x * blockDim.x;
  /*float resErr_reg[HS]; // ~ resErr_glo*/
  GPlexRegHS resErr_reg;

  for (int z = 0; z < (N-1)/grid_width  +1; z++) {
    n += z*grid_width;

    if (n < N) {

      // coordinate change
      float rotT00;
      float rotT01;
      const float r = hipo(msPar(n, 0, 0), msPar(n, 1, 0));
      rotT00 = -(msPar(n, 1, 0) + propPar(n, 1, 0))/(2*r);
      rotT01 =  (msPar(n, 0, 0) + propPar(n, 0, 0))/(2*r);

      /*float res_glo[HV];*/
      GPlexRegHV res_glo;
      subtractFirst3_fn(msPar, propPar, res_glo, N, n);

      for (int j = 0; j < HS; ++j) {
        resErr_reg[j] = 0; //resErr[j*resErr_stride + n];
      }
      addIntoUpperLeft3x3_fn(propErr, msErr, resErr_reg, N, n);

      GPlexReg2V res_loc;   //position residual in local coordinates
      RotateResidulsOnTangentPlane_fn(rotT00,rotT01,res_glo,res_loc);
      /*MPlex2S resErr_loc;//covariance sum in local position coordinates*/
      /*MPlexHH tempHH;*/
      GPlexReg2S resErr_loc; // 2x2 sym
      GPlexRegHH tempHH;  // 3*3 sym
      ProjectResErr_fn  (rotT00, rotT01, resErr_reg, tempHH);
      ProjectResErrTransp_fn(rotT00, rotT01, tempHH, resErr_loc);

      /*invertCramerSym_fn(resErr_reg);*/
      invertCramerSym2x2_fn(resErr_loc);

      chi2Similarity_fn(res_loc, resErr_loc, outChi2.ptr, outChi2.stride);
    }
  }
}


__global__ void computeChi2_kernel(
    const GPlexLS propErr, const GPlexHS msErr, const GPlexHV msPar, 
    const GPlexLV propPar, GPlexQF outChi2, const int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  if (itrack < N) {
    computeChi2_fn
      (propErr, msErr, msPar, propPar,
       outChi2, N);
  }
}


void computeChi2_wrapper(cudaStream_t &stream, 
    const GPlexLS &propErr, const GPlexHS &msErr,
    const GPlexHV &msPar, const GPlexLV &propPar, GPlexQF &outChi2,
    const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       max_blocks_x);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  computeChi2_kernel <<< grid, block, 0, stream >>>
    (propErr, msErr, msPar, propPar, outChi2, N);
 }
