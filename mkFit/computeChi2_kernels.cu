#include <stdio.h>
#include <algorithm>
#include "GPlex.h"
#include "kalmanUpdater_kernels.h"
#include "computeChi2_kernels.h"

#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

#include "HitStructuresCU.h"
#include "BinInfoUtils.h"
#include "Hit.h"

#define L 6
#define HS 6
#define HV 3
#define BLOCK_SIZE_X 32
#define MAX_BLOCKS_X 65535 // CUDA constraint


template <>
__device__ float* SVector3::ArrayCU() {
  return fArray; 
}

template <>
__device__ float* SVector6::ArrayCU() {
  return fArray; 
}

__device__ float *Hit::posArrayCU() {
  return state_.pos_.ArrayCU();
}

__device__ float *Hit::errArrayCU() {
  return state_.err_.ArrayCU();
}

__device__ void chi2Similarity_fn(
    GPlexReg2V &a,
    GPlexReg2S &c, // in registers
    float *d, size_t dN) {

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
				  float r01,//r01
				  GPlexRegHV &a  ,//res_glo
          GPlexReg2V &b  )//res_loc
{

   // res_loc = rotT * res_glo
   //   B     =  R   *    A   
  b[0] =  r00*a[0] + r01*a[1];
  b[1] =  a[2];
}

__device__ void ProjectResErr_fn(float a00,
		   float a01,
		   GPlexRegHS &b, 
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

__device__ void ProjectResErrTransp_fn(float a00,
			 float a01, GPlexRegHH &b, GPlexReg2S &c)
{
  // C = A * B, C is 3x3 sym, A is 3x3 , B is 3x3

  // Based on script generation and adapted to custom sizes.
      c[ 0] = b[ 0]*a00 + b[ 1]*a01;
      c[ 1] = b[ 3]*a00 + b[ 4]*a01;
      c[ 2] = b[ 5];
}

__device__ void computeChi2_fn(
    GPlexLS &propErr, GPlexHS &msErr, GPlexHV &msPar, GPlexLV &propPar,
    GPlexQF &outChi2, const int N) {
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
    GPlexLS propErr, GPlexHS msErr, GPlexHV msPar, GPlexLV propPar,
    GPlexQF outChi2, const int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  if (itrack < N) {
    computeChi2_fn
      (propErr, msErr, msPar, propPar,
       outChi2, N);
  }
}

void computeChi2_wrapper(cudaStream_t &stream, 
    GPlexLS &propErr, GPlexHS &msErr, // GPlex<float> resErr,
    GPlexHV &msPar, GPlexLV &propPar, GPlexQF &outChi2,
    const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  computeChi2_kernel <<< grid, block, 0, stream >>>
    (propErr, msErr, msPar, propPar, outChi2, N);
 }

template <typename GPlexObj>
__device__ void SlurpIn_fn(GPlexObj to, // float *fArray, int stride, int kSize, 
                           const char *arr, int *vi, int N) {
  int j = threadIdx.x + blockDim.x * blockIdx.x;
  if (j<N) {
    int *XHitPos = vi;
    int off = XHitPos[j] * sizeof(Hit);
    for (int i = 0; i < to.kSize; ++i) { // plex_size
      /*fArray[i*stride+ j] = * (const T*) (arr + i*sizeof(T) + off);*/
      to(j, i, 0) = * (decltype(to.ptr)) (arr + i*sizeof(decltype(*to.ptr)) + off);
    }
  }
}

template <typename GPlexObj>
__device__ void SlurpInIdx_fn(GPlexObj to, // float *fArray, int stride, int kSize, 
                             const char *arr, int idx, int N) {
  int j = threadIdx.x + blockDim.x * blockIdx.x;
  if (j<N) {
    for (int i = 0; i < to.kSize; ++i) { // plex_size
      /*fArray[i*stride+ j] = * (const T*) (arr + i*sizeof(T) + off);*/
      to(j, i, 0) = * (decltype(to.ptr)) (arr + i*sizeof(decltype(*to.ptr)) + idx);
      /*if (j == 2) {*/
        /*printf("gpu -- %d : %d, %f\n", i, idx, to(j,i,0));*/
      /*}*/
    }
  }
}


__device__ void HitToMs_fn(GPlexHS &msErr, GPlexHV &msPar,
                           Hit *hits, GPlexQI &XHitSize, GPlexHitIdx &XHitArr, 
                           int *HitsIdx, int hit_cnt, int N) {
  /*int j = threadIdx.x + blockDim.x*blockIdx.x;*/
  int itrack = threadIdx.x + blockDim.x * blockIdx.x;
  if (itrack < N) {

    const char *varr      = (char*) hits;
    const int   off_error = (char*) hits[0].errArrayCU() - varr;
    const int   off_param = (char*) hits[0].posArrayCU() - varr;

    if (hit_cnt < XHitSize[itrack]) {
      HitsIdx[itrack] = XHitArr(itrack, hit_cnt, 0) * sizeof(Hit);
    }
    SlurpInIdx_fn(msErr, varr + off_error, HitsIdx[itrack], N);
    SlurpInIdx_fn(msPar, varr + off_param, HitsIdx[itrack], N);
  }
}

__global__ void HitToMs_kernel(GPlexHS msErr, GPlexHV msPar,
    Hit *hits, GPlexQI XHitSize, GPlexHitIdx XHitArr, int *HitsIdx, int hit_cnt, int N) {

    HitToMs_fn(msErr, msPar, hits, XHitSize, XHitArr, HitsIdx, hit_cnt, N);
}

#if 1
void HitToMs_wrapper(cudaStream_t& stream,
    GPlexHS &msErr, GPlexHV &msPar, LayerOfHitsCU &layer, 
    GPlexQI &XHitSize, GPlexHitIdx &XHitArr, int *HitsIdx, int hit_cnt, int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
#if 1
  HitToMs_kernel <<< grid, block, 0 , stream >>>
    (msErr, msPar, layer.m_hits, XHitSize, XHitArr, HitsIdx, hit_cnt, N);
  cudaDeviceSynchronize();
#endif
}
#endif


__device__ void getNewBestHitChi2_fn(
    GPlexQI &XHitSize, GPlexHitIdx &XHitArr,
    float *outChi2, float &minChi2,
    int &bestHit, int hit_cnt, int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;

  if (itrack < N) {
    if (hit_cnt < XHitSize[itrack]) {
      float chi2 = fabs(outChi2[itrack]);
      if (chi2 < minChi2) {
        minChi2 = chi2;
        /*bestHit = hit_cnt;*/
        bestHit = XHitArr(itrack, hit_cnt, 0);
      }
    }
  }
}

__global__ void getNewBestHitChi2_kernel(
    GPlexQI XHitSize, GPlexHitIdx XHitArr,
    float *outChi2, float *minChi2,
    int *bestHit, int hit_cnt, int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  if (itrack < N) {
    getNewBestHitChi2_fn(XHitSize, XHitArr, outChi2, minChi2[itrack], bestHit[itrack], hit_cnt, N);
    /*printf("GPU [%d]  -- %d : %f\n", itrack, bestHit[itrack], minChi2[itrack]);*/
  }
}

void getNewBestHitChi2_wrapper(cudaStream_t &stream,
    GPlexQI &XHitSize, GPlexHitIdx &XHitArr,
    GPlexQF &outChi2, float *minChi2, int *bestHit, int hit_cnt, int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
#if 1
  getNewBestHitChi2_kernel <<< grid, block, 0, stream >>>
    (XHitSize, XHitArr, outChi2.ptr, minChi2, bestHit, hit_cnt, N);
#endif
}

void fill_array_cu(float *array, int size, float value) {
  thrust::device_ptr<float> d_ptr(array);
  thrust::fill(d_ptr, d_ptr + size, value);
}


__device__ void updateTracksWithBestHit_fn(Hit *hits, 
    float minChi2, int bestHit,
    GPlexHS &msErr, GPlexHV &msPar, GPlexLV &propPar, 
    float *Chi2, int *HitsIdx, int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  if (itrack < N) {
    /*printf("GPU [%d]  -- %d : %f\n", itrack, bestHit, minChi2);*/
    if (bestHit >= 0)
    {
      Hit   &hit  = hits[ bestHit ];
      float &chi2_local = minChi2;
	  
      /*msErr[Nhits].CopyIn(itrack, hit.errArray());*/
      /*SlurpIn_fn<float>(msErr, msErr_stride, msErr_plex_size,
        varr + (itrack*sizeof(Hit)) + off_error, XHitPos, N);*/
      /*msPar[Nhits].CopyIn(itrack, hit.posArray());*/
      /*SlurpIn_fn<float>(msPar, msPar_stride, msPar_plex_size,
        varr + (itrack*sizeof(Hit)) + off_param, XHitPos, N);*/
      for (int i = 0; i < msErr.kSize; ++i) {
        msErr(itrack, i, 0) = hit.errArrayCU()[i];
      }
      for (int i = 0; i < msPar.kSize; ++i) {
        msPar(itrack, i, 0) = hit.posArrayCU()[i];
      }
      /*Chi2(itrack, 0, 0) += chi2_local;*/
      Chi2[itrack] += chi2_local;
      /*HitsIdx[Nhits](itrack, 0, 0) = XHitPos.At(itrack, 0, 0) + bestHit[itrack];*/
      HitsIdx[itrack] = bestHit;
    }
    else
    {
      /*msErr[Nhits].SetDiagonal3x3(itrack, 666);*/
      msErr(itrack, 0, 0) = 666;
      msErr(itrack, 1, 0) = 0;
      msErr(itrack, 2, 0) = 666;
      msErr(itrack, 3, 0) = 0;
      msErr(itrack, 4, 0) = 0;
      msErr(itrack, 5, 0) = 666;

      /*msPar[Nhits](itrack,0,0) = Par[iP](itrack,0,0);*/
      /*msPar[Nhits](itrack,1,0) = Par[iP](itrack,1,0);*/
      /*msPar[Nhits](itrack,2,0) = Par[iP](itrack,2,0);*/
      for (int i = 0; i < msPar.kSize; ++i) {
        msPar(itrack, i, 0) = propPar(itrack, i, 0);
      }
      /*HitsIdx[Nhits](itrack, 0, 0) = -1;*/
      HitsIdx[itrack] = -1;

      // Don't update chi2
    }
    /*printf("GPU [%d]  -- %d : %f\n", itrack, HitsIdx[itrack], Chi2[itrack]);*/
  }
}

__global__ void updateTracksWithBestHit_kernel(Hit *hits, float *minChi2, int *bestHit,
    GPlexHS msErr, GPlexHV msPar, GPlexLV propPar, float *Chi2, int *HitsIdx, int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  if (itrack < N) {
    updateTracksWithBestHit_fn
        (hits, minChi2[itrack], bestHit[itrack],
         msErr, msPar, propPar, Chi2, HitsIdx, N);
  }
}

#if 1
void updateTracksWithBestHit_wrapper(cudaStream_t &stream,
    LayerOfHitsCU &layer, float *minChi2, int *best_hit, 
    GPlexHS &msErr, GPlexHV &msPar, GPlexLV &propPar,
    float *Chi2, int *HitsIdx, int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  updateTracksWithBestHit_kernel <<< grid, block, 0, stream >>>
      (layer.m_hits, minChi2, best_hit, msErr, msPar, propPar, Chi2, HitsIdx, N);
}
#endif

int getMaxNumHits_wrapper(GPlexQI d_XHitSize, int N) {
  thrust::device_ptr<int> d_ptr(d_XHitSize.ptr);
  int maxSize=  thrust::reduce(d_ptr, d_ptr + N, -1, thrust::maximum<int>());
  maxSize = std::min(maxSize, Config::maxHitsConsidered);

  return maxSize;
}

__global__ void bestHit_kernel(
    Hit *hits, GPlexQI XHitPos, 
    GPlexLS propErr, GPlexHS msErr, GPlexHV msPar,
    GPlexLV propPar, GPlexQF outChi2,
    /*float* propErr, size_t propErr_stride,*/
    /*float* msErr, size_t msErr_stride, size_t msErr_plex_size,*/
    /*float *msPar, size_t msPar_stride, size_t msPar_plex_size,*/
    /*float *propPar, size_t propPar_stride,*/
    /*float *outChi2, size_t outChi2_stride,*/
    float *Chi2, int *HitsIdx,
    int maxSize, int N) {

  /*int itrack = threadIdx.x + blockDim.x*blockIdx.x;*/
  int bestHit_reg = -1;
  float minChi2_reg = 15.f;

  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
    /*HitToMs_fn(msErr, msPar, hits, XHitPos, hit_cnt, N);*/
#if 0
      // TODO: add CMSGeom
      if (Config::useCMSGeom) {
        //propagateHelixToRMPlex(psErr,  psPar, inChg,  msPar, propErr, propPar);
        throw std::runtime_error("useCMSGeom not implemented yet for GPU");
      } else {}
#endif
    computeChi2_fn(propErr, msErr, msPar, propPar, outChi2, N);
    /*getNewBestHitChi2_fn(outChi2.ptr, minChi2_reg, bestHit_reg, hit_cnt, N);*/
  }
  /*updateTracksWithBestHit_fn*/
      /*(hits, XHitPos,*/
       /*minChi2_reg, bestHit_reg,*/
       /*msErr, msPar, propPar,*/
       /*Chi2, HitsIdx,*/
       /*N);*/
}

#if 0
void bestHit_wrapper(cudaStream_t &stream,
    BunchOfHitsCU &bunch, GPlexQI &XHitPos, 
    GPlexLS &propErr, GPlexHS &msErr, GPlexHV &msPar,
    GPlexLV &propPar, GPlexQF &outChi2,
    float *Chi2, int *HitsIdx,
    int maxSize, int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  bestHit_kernel <<< grid, block, 0, stream >>>
    (bunch.m_hits, XHitPos,
     propErr, msErr, msPar, propPar, outChi2,
     /*propErr.ptr, propErr.stride,*/
     /*msErr.ptr, msErr.stride, msErr.kSize,*/
     /*msPar.ptr, msPar.stride, msPar.kSize,*/
     /*outChi2.ptr, outChi2.stride,*/
     Chi2, HitsIdx,
     maxSize, N);
}
#endif

__global__ void selectHitRanges_kernel(Hit *hits,
    int *phi_bin_infos_first, int *phi_bin_infos_second, int bunch_fill_index,
    GPlexQI XHitPos, GPlexQI XHitSize, GPlexLS Err, GPlexLV Par,
    bool useCMSGeom, int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  if (itrack < N) {
    // must store hit vector into a data member so it can be used in hit selection.
    // or ... can have it passed into other functions.
    // somewhat yucky, either way.

    // Also, must store two ints per Matriplex elements ... first index and size.
    // These are XPos and XSize

    /*const int iI = iP;*/
    // Hmmh ... this should all be solved by partitioning ... let's try below ...
    //
    // float eta = getEta(eta_predx,eta_predy,eta_predz);
    // //protect against anomalous eta (should go into getEtaPartition maybe?)
    // if (fabs(eta) > etaDet) eta = (eta>0 ? etaDet*0.99 : -etaDet*0.99);
    // unsigned int etabin = getEtaPartition(eta,etaDet);

    const float predx = Par(itrack, (0*1 + 0), 0);  // Par[iI].ConstAt(itrack, 0, 0);
    const float predy = Par(itrack, (1*1 + 0), 0);  // Par[iI].ConstAt(itrack, 1, 0);
    const float predz = Par(itrack, (2*1 + 0), 0);  // Par[iI].ConstAt(itrack, 2, 0);

    float phi = getPhi(predx,predy);

    const float px2py2 = predx*predx+predy*predy; // predicted radius^2
    const float dphidx = -predy/px2py2;
    const float dphidy =  predx/px2py2;
    // const float dphi2  =     dphidx*dphidx*(Err[iI].ConstAt(itrack, 0, 0) /*propState.errors.At(0,0)*/) +
    //                          dphidy*dphidy*(Err[iI].ConstAt(itrack, 1, 1) /*propState.errors.At(1,1)*/) +
    //                      2 * dphidx*dphidy*(Err[iI].ConstAt(itrack, 0, 1) /*propState.errors.At(0,1)*/);
    const float dphi2  =     dphidx*dphidx*Err(itrack, 0, 0) +
                             dphidy*dphidy*Err(itrack, 2, 0) +
                         2 * dphidx*dphidy*Err(itrack, 1, 0);

    const float dphi       = sqrtf(fabs(dphi2));//how come I get negative squared errors sometimes? MT -- how small?
    const float nSigmaDphi = fminf(fmaxf(Config::nSigma*dphi, Config::minDPhi), Config::PI);
    //const float nSigmaDphi = Config::nSigma*dphi;

    float dPhiMargin = 0.;
    if (useCMSGeom) {
      //now correct for bending and for layer thickness unsing linear approximation
      /*const float predpx = Par[iP].ConstAt(itrack, 3, 0);*/
      /*const float predpy = Par[iP].ConstAt(itrack, 4, 0);*/
      const float predpx = Par(itrack, (3*1 + 0), 0);
      const float predpy = Par(itrack, (4*1 + 0), 0);
      float deltaR = Config::cmsDeltaRad; //fixme! using constant vale, to be taken from layer properties
      float radius = sqrt(px2py2);
      float pt     = sqrt(predpx*predpx + predpy*predpy);
      float cosTheta = ( predx*predpx + predy*predpy )/(pt*radius);
      float hipo = deltaR/cosTheta;
      float dist = sqrt(hipo*hipo - deltaR*deltaR);
      dPhiMargin = dist/radius;
    }
    const float dphiMinus = normalizedPhi(phi-nSigmaDphi-dPhiMargin);
    const float dphiPlus  = normalizedPhi(phi+nSigmaDphi+dPhiMargin);
// FIXME ^ OK

#ifdef DEBUG
    std::ostringstream xout;
    bool               xout_dump = false;
    xout << "--------------------------------------------------------------------------------\n";
    xout << "phi  = " << phi  << ", dphiMinus = " << dphiMinus << ", dphiPlus = " << dphiPlus << std::endl;
    xout << "dphi = " << dphi  << ", dphi2 = " << dphi2 << ", nSigmaDphi = " << nSigmaDphi << ", nSigma = " << Config::nSigma << std::endl;
#endif

    int   phiBinMinus = getPhiPartition(dphiMinus);
    int   phiBinPlus  = getPhiPartition(dphiPlus);

#ifdef DEBUG
    xout << "phiBinMinus = " << phiBinMinus << ", phiBinPlus = " << phiBinPlus << std::endl;
#endif

    // XXXX are these checks really needed?
    phiBinMinus = fmaxf(0,phiBinMinus);
    phiBinMinus = fminf(Config::nPhiPart-1,phiBinMinus);
    phiBinPlus  = fmaxf(0,phiBinPlus);
    phiBinPlus  = fminf(Config::nPhiPart-1,phiBinPlus);

    //PhiBinInfo_t binInfoMinus = bunch_of_hits.m_phi_bin_infos[phiBinMinus];
    //PhiBinInfo_t binInfoPlus  = bunch_of_hits.m_phi_bin_infos[phiBinPlus];
    int binInfoMinus_first = phi_bin_infos_first[phiBinMinus];
    int binInfoMinus_second = phi_bin_infos_second[phiBinMinus];
    int binInfoPlus_first = phi_bin_infos_first[phiBinPlus];
    int binInfoPlus_second = phi_bin_infos_second[phiBinPlus];


    /*if (binInfoPlus.first + binInfoPlus.second - binInfoMinus.first > Config::maxHitsConsidered)*/
    if (binInfoPlus_first + binInfoPlus_second - binInfoMinus_first > Config::maxHitsConsidered)
    {
      // XXXX
      // Do something smart to reduce the range.
      // I'd go for taking the exact phi bin and then walking left and right ...
      // but this gives the wrap-around problem again.
    }

    // XXXX
    // Hmmh ... maybe the copying of extras should be done on demand.
    // BunchOfHits could know how many extras it has already.
    // Or Giuseppe is right ... and we should just update the index vector for SlurpIn
    // instead of shifting of the base address as is done now. Sigh.
    
    // fixme: temporary to avoid wrapping
    // This is now fixed with Config::maxHitsConsidered extra hits copied to the end +
    // changing XHitBegin/End to XHitPos/Size.
    // Putting all of it into DEBUG
#ifdef DEBUG
    if (binInfoMinus > binInfoPlus)
    {
      // xout_dump = true;
      xout << "FIXER IN:  phiBinMinus = " << phiBinMinus << ", phiBinPlus = " << phiBinPlus << std::endl;
      xout << "FIXER IN:  BIMinus.first = " << binInfoMinus.first << ", BIPlus.first = " << binInfoPlus.first << std::endl;
      xout << "FIXER IN:  BIMinus.second = " << binInfoMinus.second << ", BIPlus.second = " << binInfoPlus.second << std::endl;

      int phibin = getPhiPartition(phi);

      xout << "FIXER   :  phibin = " << phibin << std::endl;

      // XXXX are those two really needed?
      phibin = std::max(0,phibin);
      phibin = std::min(Config::nPhiPart-1,phibin);

      xout << "FIXER   :  phibin = " << phibin << std::endl;
    }
#endif

    XHitPos[itrack] = binInfoMinus_first;
    XHitSize[itrack] = binInfoPlus_first + binInfoPlus_second - binInfoMinus_first;
    if (XHitSize[itrack] < 0)
    {
      // XXX It would be nice to have BunchOfHits.m_n_real_hits.
      /*XHitSize[itrack] += bunch_of_hits.m_fill_index - Config::maxHitsConsidered;*/
      XHitSize[itrack] += bunch_fill_index - Config::maxHitsConsidered;
    }

    // XXXX Hack to limit N_hits to maxHitsConsidered.
    // Should at least take hits around central bin -- to be explored, esp. with jet presence.
    // Strange ... this is worse than just taking first 25 hits !!!
    // Comment out for now. Must talk to Giuseppe about this.
    // if (XHitSize.At(itrack, 0, 0) > Config::maxHitsConsidered)
    // {
    //   xout_dump = true;
    //   XHitPos .At(itrack, 0, 0) += (XHitSize.At(itrack, 0, 0) - Config::maxHitsConsidered) / 2;
    //   XHitSize.At(itrack, 0, 0) = Config::maxHitsConsidered;
    // }

#ifdef DEBUG
    xout << "found range firstHit=" << XHitPos.At(itrack, 0, 0) << " size=" << XHitSize.At(itrack, 0, 0) << std::endl;
    if (xout_dump)
       std::cout << xout.str();
#endif

  }
}

#if 0
void selectHitRanges_wrapper(cudaStream_t &stream, BunchOfHitsCU &bunch, 
    GPlexQI &XHitPos, GPlexQI &XHitSize,
    GPlexLS &Err, GPlexLV &Par,
    int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  selectHitRanges_kernel <<< grid, block, 0, stream >>>
    (bunch.m_hits, bunch.m_phi_bin_infos_first, 
     bunch.m_phi_bin_infos_second, bunch.m_fill_index,
     XHitPos, XHitSize, Err, Par,
     Config::useCMSGeom, N);
}
#endif
