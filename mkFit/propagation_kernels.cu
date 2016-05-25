#include "Config.h"
#include "Debug.h"
#include "propagation_kernels.h"
#include <stdio.h>

constexpr int L = 6;
constexpr int LL2 = 36;
constexpr int LS = 21;

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

// values from 32 to 512 give good results.
// 32 gives slightly better results (on a K40)
#define BLOCK_SIZE_X 32
#define MAX_BLOCKS_X 65535 // CUDA constraint

#if 0
__device__ float hipo(float x, float y) {
  return std::sqrt(x*x + y*y);
}

__device__ void sincos4_cu(float x, float& sin, float& cos) {
   // Had this writen with explicit division by factorial.
   // The *whole* fitting test ran like 2.5% slower on MIC, sigh.
   cos  = 1;
   sin  = x;   x *= x * 0.5f;
   cos -= x;   x *= x * 0.33333333f;
   sin -= x;   x *= x * 0.25f;
   cos += x;
}
#endif

// computeJacobianSimple works on values that are in registers.
// Registers are thread-private. Thus this function has no notion of
// parallelism. It is ran serially by each calling thread.
__device__ void computeJacobianSimple(float *errorProp,
    float s, float k, float p, float pxin, float pyin, float pzin, 
    float TP, float cosTP, float sinTP, int N) {

  // std::cout << "total path s=" << s << std::endl;
  // TD = s*pt/p;
  // TP = TD/(pt*k) = s/(p*k);
  float dTPdpx = -s*pxin/(k*p*p*p);
  float dTPdpy = -s*pyin/(k*p*p*p);
  float dTPdpz = -s*pzin/(k*p*p*p);
  //ok let's assume that the quantity with no error is the angular path (phase change)
  dTPdpx = 0;
  dTPdpy = 0;
  dTPdpz = 0;
  
  //derive these to compute jacobian
  //x = xin + k*(pxin*sinTP-pyin*(1-cosTP));
  //y = yin + k*(pyin*sinTP+pxin*(1-cosTP));
  //z = zin + k*TP*pzin;
  //px = pxin*cosTP-pyin*sinTP;
  //py = pyin*cosTP+pxin*sinTP;
  //pz = pzin;
  //jacobian
  
  errorProp[(0*L + 0)] = 1.;	                                             //dxdx
  errorProp[(0*L + 1)] = 0.;	                                             //dxdy
  errorProp[(0*L + 2)] = 0.;                                                     //dxdz
  errorProp[(0*L + 3)] = k*(sinTP + pxin*cosTP*dTPdpx - pyin*sinTP*dTPdpx);      //dxdpx
  errorProp[(0*L + 4)] = k*(pxin*cosTP*dTPdpy - 1. + cosTP - pyin*sinTP*dTPdpy); //dxdpy
  errorProp[(0*L + 5)] = k*dTPdpz*(pxin*cosTP - pyin*sinTP);                     //dxdpz
  errorProp[(1*L + 0)] = 0.;	                                             //dydx
  errorProp[(1*L + 1)] = 1.;	                                             //dydy
  errorProp[(1*L + 2)] = 0.;                                                     //dydz
  errorProp[(1*L + 3)] = k*(pyin*cosTP*dTPdpx + 1. - cosTP + pxin*sinTP*dTPdpx); //dydpx
  errorProp[(1*L + 4)] = k*(sinTP + pyin*cosTP*dTPdpy + pxin*sinTP*dTPdpy);      //dydpy
  errorProp[(1*L + 5)] = k*dTPdpz*(pyin*cosTP + pxin*sinTP);                     //dydpz
  errorProp[(2*L + 0)] = 0.;	                                             //dzdx
  errorProp[(2*L + 1)] = 0.;	                                             //dzdy
  errorProp[(2*L + 2)] = 1.;                                                     //dzdz
  errorProp[(2*L + 3)] = k*pzin*dTPdpx;                                          //dzdpx
  errorProp[(2*L + 4)] = k*pzin*dTPdpy;                                          //dzdpy
  errorProp[(2*L + 5)] = k*(TP + dTPdpz*pzin);                                   //dzdpz
  errorProp[(3*L + 0)] = 0.;	                                             //dpxdx
  errorProp[(3*L + 1)] = 0.;	                                             //dpxdy
  errorProp[(3*L + 2)] = 0.;                                                     //dpxdz
  errorProp[(3*L + 3)] = cosTP - dTPdpx*(pxin*sinTP + pyin*cosTP);               //dpxdpx
  errorProp[(3*L + 4)] = -sinTP - dTPdpy*(pxin*sinTP + pyin*cosTP);              //dpxdpy
  errorProp[(3*L + 5)] = -dTPdpz*(pxin*sinTP + pyin*cosTP);                      //dpxdpz
  errorProp[(4*L + 0)] = 0.;                                                     //dpydx
  errorProp[(4*L + 1)] = 0.;	                                             //dpydy
  errorProp[(4*L + 2)] = 0.;                                                     //dpydz
  errorProp[(4*L + 3)] = +sinTP - dTPdpx*(pyin*sinTP - pxin*cosTP);              //dpydpx
  errorProp[(4*L + 4)] = +cosTP - dTPdpy*(pyin*sinTP - pxin*cosTP);              //dpydpy
  errorProp[(4*L + 5)] = -dTPdpz*(pyin*sinTP - pxin*cosTP);                      //dpydpz
  errorProp[(5*L + 0)] = 0.;                                                     //dpzdx
  errorProp[(5*L + 1)] = 0.;						     //dpzdy
  errorProp[(5*L + 2)] = 0.;						     //dpzdz 
  errorProp[(5*L + 3)] = 0.;						     //dpzdpx
  errorProp[(5*L + 4)] = 0.;						     //dpzdpy
  errorProp[(5*L + 5)] = 1.;						     //dpzdpz  
}

/// Compute MsRad /////////////////////////////////////////////////////////////
__device__ void assignMsRad_fn(const float r, float* msRad, int N, int n) {
  /*int n = threadIdx.x + blockIdx.x * blockDim.x;*/
  if (n < N) {
    *msRad = r;
  }
}

// Not passing msRad.stride, as QF == 1 (second dim f msRad)
__device__ void computeMsRad_fn(const GPlexHV& __restrict__ msPar,
    float* msRad, int N, int n) {
  /*int n = threadIdx.x + blockIdx.x * blockDim.x;*/
  if (n < N) {
    *msRad = hipo(msPar.ptr[n], msPar.ptr[n + msPar.stride]);
  }
}

#include "PropagationMPlex.icc"

__device__ 
void helixAtRFromIterative_fn(const GPlexLV& inPar,
    const GPlexQI& inChg, GPlexLV& outPar_global, const GPlexReg<float,1,1>& msRad, 
    GPlexReg<float, LL2, L>& errorProp, int N, int n) {

  /*int n = threadIdx.x + blockIdx.x * blockDim.x;*/

  GPlexReg<float, LL2, 1> outPar;

  if (n < N) {
    for (int j = 0; j < 5; ++j) {
      outPar[j] = outPar_global(n, j, 0);
    }
    errorProp.SetVal(0);

    helixAtRFromIterative_impl(inPar, inChg, outPar, msRad, errorProp, n, n+1);

    // Once computations are done. Get values from registers to global memory.
    for (int j = 0; j < 5; ++j) {
      outPar_global(n, j, 0) = outPar[j];
    }
  }
}

/// Similarity ////////////////////////////////////////////////////////////////
__device__ void similarity_fn(float* a, float *b, size_t stride_outErr,
    int N, int n) {
  size_t bN = stride_outErr;
  
  // Keep most values in registers.
  float b_reg[LL2];
  // To avoid using too many registers, tmp[] as a limited size and is reused.
  float tmp[6];

  /*int n = threadIdx.x + blockIdx.x * blockDim.x;*/

  if (n < N) {
    for (int j = 0; j < LS; j++) {
      b_reg[j] = b[n + j*bN];
    }

    tmp[ 0] = a[0]*b_reg[ 0] + a[1]*b_reg[ 1] + a[3]*b_reg[ 6] + a[4]*b_reg[10];
    tmp[ 1] = a[0]*b_reg[ 1] + a[1]*b_reg[ 2] + a[3]*b_reg[ 7] + a[4]*b_reg[11];
    /*tmp[ 2] = a[0]*b_reg[ 3] + a[1]*b_reg[ 4] + a[3]*b_reg[ 8] + a[4]*b_reg[12];*/
    tmp[ 3] = a[0]*b_reg[ 6] + a[1]*b_reg[ 7] + a[3]*b_reg[ 9] + a[4]*b_reg[13];
    tmp[ 4] = a[0]*b_reg[10] + a[1]*b_reg[11] + a[3]*b_reg[13] + a[4]*b_reg[14];
    /*tmp[ 5] = a[0]*b_reg[15] + a[1]*b_reg[16] + a[3]*b_reg[18] + a[4]*b_reg[19];*/

    b[ 0*bN+n] = tmp[ 0]*a[0] + tmp[ 1]*a[1] + tmp[ 3]*a[3] + tmp[ 4]*a[4];


    tmp[ 0] = a[6]*b_reg[ 0] + a[7]*b_reg[ 1] + a[9]*b_reg[ 6] + a[10]*b_reg[10];
    tmp[ 1] = a[6]*b_reg[ 1] + a[7]*b_reg[ 2] + a[9]*b_reg[ 7] + a[10]*b_reg[11];
    /*tmp[ 8] = a[6]*b_reg[ 3] + a[7]*b_reg[ 4] + a[9]*b_reg[ 8] + a[10]*b_reg[12];*/
    tmp[ 3] = a[6]*b_reg[ 6] + a[7]*b_reg[ 7] + a[9]*b_reg[ 9] + a[10]*b_reg[13];
    tmp[ 4] = a[6]*b_reg[10] + a[7]*b_reg[11] + a[9]*b_reg[13] + a[10]*b_reg[14];
    /*tmp[11] = a[6]*b_reg[15] + a[7]*b_reg[16] + a[9]*b_reg[18] + a[10]*b_reg[19];*/

    b[ 1*bN+n] = tmp[ 0]*a[0] + tmp[ 1]*a[1] + tmp[ 3]*a[3] + tmp[ 4]*a[4];
    b[ 2*bN+n] = tmp[ 0]*a[6] + tmp[ 1]*a[7] + tmp[ 3]*a[9] + tmp[ 4]*a[10];


    tmp[ 0] = a[12]*b_reg[ 0] + a[13]*b_reg[ 1] + b_reg[ 3] + a[15]*b_reg[ 6] + a[16]*b_reg[10] + a[17]*b_reg[15];
    tmp[ 1] = a[12]*b_reg[ 1] + a[13]*b_reg[ 2] + b_reg[ 4] + a[15]*b_reg[ 7] + a[16]*b_reg[11] + a[17]*b_reg[16];
    tmp[ 2] = a[12]*b_reg[ 3] + a[13]*b_reg[ 4] + b_reg[ 5] + a[15]*b_reg[ 8] + a[16]*b_reg[12] + a[17]*b_reg[17];
    tmp[ 3] = a[12]*b_reg[ 6] + a[13]*b_reg[ 7] + b_reg[ 8] + a[15]*b_reg[ 9] + a[16]*b_reg[13] + a[17]*b_reg[18];
    tmp[ 4] = a[12]*b_reg[10] + a[13]*b_reg[11] + b_reg[12] + a[15]*b_reg[13] + a[16]*b_reg[14] + a[17]*b_reg[19];
    tmp[ 5] = a[12]*b_reg[15] + a[13]*b_reg[16] + b_reg[17] + a[15]*b_reg[18] + a[16]*b_reg[19] + a[17]*b_reg[20];

    b[ 3*bN+n] = tmp[ 0]*a[0] + tmp[ 1]*a[1]           + tmp[ 3]*a[3] + tmp[ 4]*a[4];
    b[ 4*bN+n] = tmp[ 0]*a[6] + tmp[ 1]*a[7]           + tmp[ 3]*a[9] + tmp[ 4]*a[10];
    b[ 5*bN+n] = tmp[ 0]*a[12] + tmp[ 1]*a[13] + tmp[ 2] + tmp[ 3]*a[15] + tmp[ 4]*a[16] + tmp[ 5]*a[17];


    tmp[ 0] = a[18]*b_reg[ 0] + a[19]*b_reg[ 1] + a[21]*b_reg[ 6] + a[22]*b_reg[10];
    tmp[ 1] = a[18]*b_reg[ 1] + a[19]*b_reg[ 2] + a[21]*b_reg[ 7] + a[22]*b_reg[11];
    tmp[ 2] = a[18]*b_reg[ 3] + a[19]*b_reg[ 4] + a[21]*b_reg[ 8] + a[22]*b_reg[12];
    tmp[ 3] = a[18]*b_reg[ 6] + a[19]*b_reg[ 7] + a[21]*b_reg[ 9] + a[22]*b_reg[13];
    tmp[ 4] = a[18]*b_reg[10] + a[19]*b_reg[11] + a[21]*b_reg[13] + a[22]*b_reg[14];
    tmp[ 5] = a[18]*b_reg[15] + a[19]*b_reg[16] + a[21]*b_reg[18] + a[22]*b_reg[19];

    b[ 6*bN+n] = tmp[ 0]*a[0] + tmp[ 1]*a[1]           + tmp[ 3]*a[3] + tmp[ 4]*a[4];
    b[ 7*bN+n] = tmp[ 0]*a[6] + tmp[ 1]*a[7]           + tmp[ 3]*a[9] + tmp[ 4]*a[10];
    b[ 8*bN+n] = tmp[ 0]*a[12] + tmp[ 1]*a[13] + tmp[ 2] + tmp[ 3]*a[15] + tmp[ 4]*a[16] + tmp[ 5]*a[17];
    b[ 9*bN+n] = tmp[ 0]*a[18] + tmp[ 1]*a[19]           + tmp[ 3]*a[21] + tmp[ 4]*a[22];


    tmp[ 0] = a[24]*b_reg[ 0] + a[25]*b_reg[ 1] + a[27]*b_reg[ 6] + a[28]*b_reg[10];
    tmp[ 1] = a[24]*b_reg[ 1] + a[25]*b_reg[ 2] + a[27]*b_reg[ 7] + a[28]*b_reg[11];
    tmp[ 2] = a[24]*b_reg[ 3] + a[25]*b_reg[ 4] + a[27]*b_reg[ 8] + a[28]*b_reg[12];
    tmp[ 3] = a[24]*b_reg[ 6] + a[25]*b_reg[ 7] + a[27]*b_reg[ 9] + a[28]*b_reg[13];
    tmp[ 4] = a[24]*b_reg[10] + a[25]*b_reg[11] + a[27]*b_reg[13] + a[28]*b_reg[14];
    tmp[ 5] = a[24]*b_reg[15] + a[25]*b_reg[16] + a[27]*b_reg[18] + a[28]*b_reg[19];

    b[10*bN+n] = tmp[ 0]*a[0] + tmp[ 1]*a[1]           + tmp[ 3]*a[3] + tmp[ 4]*a[4];
    b[11*bN+n] = tmp[ 0]*a[6] + tmp[ 1]*a[7]           + tmp[ 3]*a[9] + tmp[ 4]*a[10];
    b[12*bN+n] = tmp[ 0]*a[12] + tmp[ 1]*a[13] + tmp[ 2] + tmp[ 3]*a[15] + tmp[ 4]*a[16] + tmp[ 5]*a[17];
    b[13*bN+n] = tmp[ 0]*a[18] + tmp[ 1]*a[19]           + tmp[ 3]*a[21] + tmp[ 4]*a[22];
    b[14*bN+n] = tmp[ 0]*a[24] + tmp[ 1]*a[25]           + tmp[ 3]*a[27] + tmp[ 4]*a[28];

    tmp[ 0] = b_reg[15];
    tmp[ 1] = b_reg[16];
    tmp[ 2] = b_reg[17];
    tmp[ 3] = b_reg[18];
    tmp[ 4] = b_reg[19];
    tmp[ 5] = b_reg[20];

    // MultHelixPropTransp
    b[15*bN+n] = tmp[ 0]*a[0] + tmp[ 1]*a[1]           + tmp[ 3]*a[3] + tmp[ 4]*a[4];
    b[16*bN+n] = tmp[ 0]*a[6] + tmp[ 1]*a[7]           + tmp[ 3]*a[9] + tmp[ 4]*a[10];
    b[17*bN+n] = tmp[ 0]*a[12] + tmp[ 1]*a[13] + tmp[ 2] + tmp[ 3]*a[15] + tmp[ 4]*a[16] + tmp[ 5]*a[17];
    b[18*bN+n] = tmp[ 0]*a[18] + tmp[ 1]*a[19]           + tmp[ 3]*a[21] + tmp[ 4]*a[22];
    b[19*bN+n] = tmp[ 0]*a[24] + tmp[ 1]*a[25]           + tmp[ 3]*a[27] + tmp[ 4]*a[28];
    b[20*bN+n] = tmp[ 5];
  }
}

// PropagationMPlex.cc:propagateHelixToRMPlex, first version with 6 arguments 
__global__ void propagation_kernel(
    GPlexHV msPar,
    GPlexLV inPar, GPlexQI inChg,
    GPlexLV outPar, GPlexLL errorProp,
    GPlexLS outErr, int N) {

  int grid_width = blockDim.x * gridDim.x;
  int n = threadIdx.x + blockIdx.x * blockDim.x;
  GPlexReg<float,1,1> msRad_reg;
  // Using registers instead of shared memory is ~ 30% faster.
  GPlexReg<float, LL2, L> errorProp_reg;
  // If there is more matrices than MAX_BLOCKS_X * BLOCK_SIZE_X 
  for (int z = 0; z < (N-1)/grid_width  +1; z++) {
    n += z*grid_width;
    if (n < N) {
#if 0
      computeMsRad_fn(msPar, stride_msPar, &msRad_reg, N, n);
      if (Config::doIterative) {
        helixAtRFromIterative_fn(inPar, inPar_stride,
            inChg, outPar, outPar_stride, msRad_reg, 
            errorProp_reg, N, n);
      } else {
        // TODO: not ported for now. Assuming Config::doIterative
        // helixAtRFromIntersection(inPar, inChg, outPar, msRad, errorProp);
      }
      similarity_fn(errorProp_reg, outErr, outErr_stride, N, n);
#endif
      computeMsRad_fn(msPar, msRad_reg.arr, N, n);
      helixAtRFromIterative_fn(inPar, inChg, outPar, msRad_reg, errorProp_reg, N, n);
      similarity_fn(errorProp_reg.arr, outErr.ptr, outErr.stride, N, n);
    }
  }
}


void propagation_wrapper(cudaStream_t& stream,
    GPlexHV& msPar,
    GPlexLV& inPar, GPlexQI& inChg,
    GPlexLV& outPar, GPlexLL& errorProp,
    GPlexLS& outErr, 
    const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  propagation_kernel <<<grid, block, 0, stream >>>(msPar, inPar, inChg, outPar, errorProp, outErr, N);
}


// PropagationMPlex.cc:propagateHelixToRMPlex, second version with 7 arguments 
// Imposes the radius
__global__ void propagationForBuilding_kernel(
    float r,
    float *inPar, size_t inPar_stride, int *inChg,
    float *outPar, size_t outPar_stride, float *errorProp,
    size_t errorProp_stride, float *outErr, size_t outErr_stride, int N) {
#if 0
  int grid_width = blockDim.x * gridDim.x;
  int n = threadIdx.x + blockIdx.x * blockDim.x;
  float msRad_reg;
  // Using registers instead of shared memory is ~ 30% faster.
  float errorProp_reg[LL2];
  // If there is more matrices than MAX_BLOCKS_X * BLOCK_SIZE_X 
  for (int z = 0; z < (N-1)/grid_width  +1; z++) {
    n += z*grid_width;
    if (n < N) {
      assignMsRad_fn(r, &msRad_reg, N, n);
      if (Config::doIterative) {
        helixAtRFromIterative_fn(inPar, inPar_stride,
            inChg, outPar, outPar_stride, msRad_reg, 
            errorProp_reg, N, n);
      } else {
        // TODO: not ported for now. Assuming Config::doIterative
        // helixAtRFromIntersection(inPar, inChg, outPar, msRad, errorProp);
      }
      similarity_fn(errorProp_reg, outErr, outErr_stride, N, n);
    }
  }
#endif
}

void propagationForBuilding_wrapper(cudaStream_t& stream,
    float radius,
    GPlexLV& inPar, GPlexQI& inChg,
    GPlexLV& outPar, GPlexLL& errorProp,
    GPlexLS& outErr, 
    const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  propagationForBuilding_kernel<<<grid, block, 0, stream >>>
    (radius,
     inPar.ptr, inPar.stride, inChg.ptr,
     outPar.ptr, outPar.stride, errorProp.ptr,
     errorProp.stride, outErr.ptr, outErr.stride, N);
}

