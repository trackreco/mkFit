#include <stdio.h>
#include <algorithm>
#include "GPlex.h"
#include "kalmanUpdater_kernels.h"
#include "computeChi2_kernels.h"
#include "HitStructuresCU.h"

#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

#include "Hit.h"

#define L 6
#define HS 6
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
    float *a, size_t aN,
    float *b, size_t bN,
    float *c, // in registers
    float *d, size_t dN) {

  int n = threadIdx.x + blockIdx.x * blockDim.x;

  // manually subrtact into local vars -- 3 of them
  float x0 = a[0 * aN + n] - b[0 * aN + n];
  float x1 = a[1 * aN + n] - b[1 * aN + n];
  float x2 = a[2 * aN + n] - b[2 * aN + n];
  d[0 * dN + n] = c[0]*x0*x0 + c[2]*x1*x1 + c[5]*x2*x2 +
              2*( c[1]*x1*x0 + c[3]*x2*x0 + c[4]*x1*x2);
}

__device__ void computeChi2_fn(
    float* propErr, size_t propErr_stride,
    float* msErr, size_t msErr_stride,
    /*float* resErr, size_t resErr_stride,*/
    float *msPar, size_t msPar_stride,
    float *propPar, size_t propPar_stride,
    float *outChi2, size_t outChi2_stride,
    const int N) {
  int grid_width = blockDim.x * gridDim.x;
  int n = threadIdx.x + blockIdx.x * blockDim.x;
  float resErr_reg[HS];

  for (int z = 0; z < (N-1)/grid_width  +1; z++) {
    n += z*grid_width;
    if (n < N) {
      for (int j = 0; j < HS; ++j) {
        resErr_reg[j] = 0; //resErr[j*resErr_stride + n];
      }
      addIntoUpperLeft3x3_fn(propErr, propErr_stride,
          msErr, msErr_stride, resErr_reg, N, n);
      invertCramerSym_fn(resErr_reg);

      chi2Similarity_fn(msPar, msPar_stride,
          propPar, propPar_stride, resErr_reg, 
          outChi2, outChi2_stride);
      /*for (int j = 0; j < HS; ++j) {*/
        /*resErr[j*resErr_stride + n] = resErr_reg[j];*/
      /*}*/
    }
  }
}

void computeChi2_wrapper(cudaStream_t &stream, 
    GPlexLS propErr, GPlexHS msErr, // GPlex<float> resErr,
    GPlexHV msPar, GPlexLV propPar, GPlexQF outChi2,
    const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
#if 0
  computeChi2_kernel <<< grid, block, 0, stream >>>
    (propErr.ptr, propErr.stride,
     msErr.ptr, msErr.stride,
     /*resErr.ptr, resErr.stride,*/
     msPar.ptr, msPar.stride,
     propPar.ptr, propPar.stride,
     outChi2.ptr, outChi2.stride,
     N);
#endif
 }

template <typename T>
__device__ void SlurpIn_fn(float *fArray, int stride, int kSize, 
                           const char *arr, int *vi, int N) {
  int j = threadIdx.x + blockDim.x * blockIdx.x;
  if (j<N) {
    for (int i = 0; i < kSize; ++i) { // plex_size
      int *XHitPos = vi;
      int off = XHitPos[j] * sizeof(Hit);
      fArray[i*stride+ j] = * (const T*) (arr + i*sizeof(T) + off);
    }
  }
}


__device__ void HitToMs_fn(float *msErr, int msErr_stride, int msErr_plex_size,
                           float *msPar, int msPar_stride, int msPar_plex_size,
                           Hit *hits, int *XHitPos, int hit_cnt, int N) {
  /*int j = threadIdx.x + blockDim.x*blockIdx.x;*/

  const char *varr      = (char*) hits;
  const int   off_error = (char*) hits[0].errArrayCU() - varr;
  const int   off_param = (char*) hits[0].posArrayCU() - varr;
  
  /*if (j<N) {*/
  /*if (j==1) {*/
    /*int hi = XHitPos[j];*/
    /*Hit *hits_shifted = &hits[hit_cnt];*/
    /*for (int i = 0; i < msPar_plex_size; ++i) {*/
      /*msPar[i*msPar_stride+ j] = hits_shifted[hi].posArrayCU()[i];*/
      /*if (msPar[i*msPar_stride+ j] != hits_shifted[hi].posArrayCU()[i]) {*/
        /*printf("(C:%f/G:%f)    ", msPar[i*msPar_stride + j], hits_shifted[hi].posArrayCU()[i]);*/
      /*}*/
    /*}*/
    /*printf("\n");*/
  /*}*/

  SlurpIn_fn<float>(msErr, msErr_stride, msErr_plex_size, varr + (hit_cnt*sizeof(Hit)) + off_error, XHitPos, N);
  SlurpIn_fn<float>(msPar, msPar_stride, msPar_plex_size, varr + (hit_cnt*sizeof(Hit)) + off_param, XHitPos, N);
  
  /*if (j==2) {*/
    /*for (int i = 0; i < msPar_plex_size; ++i) {*/
      /*printf("GPU:Par[%d*N+%d]=%f   ", i, j, msPar[i*msErr_stride + j]);*/
    /*}*/
  /*}*/
}


void HitToMs_wrapper(cudaStream_t& stream,
    GPlexHS &msErr, GPlexHV &msPar, BunchOfHitsCU &bunch, 
    GPlexQI &XHitPos, int hit_cnt, int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
#if 0
  HitToMs_kernel <<< grid, block, 0 , stream >>>
    (msErr.ptr, msErr.stride, msErr.y,
     msPar.ptr, msPar.stride, msPar.y,
     bunch.m_hits, XHitPos.ptr, hit_cnt, N);
#endif
}


__device__ void getNewBestHitChi2_fn(float *outChi2, float &minChi2,
    int &bestHit, int hit_cnt, int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;

  if (itrack < N) {
    float chi2 = fabs(outChi2[itrack]);
    if (chi2 < minChi2) {
      minChi2 = chi2;
      bestHit = hit_cnt;
    }
  }
}

void getNewBestHitChi2_wrapper(cudaStream_t &stream,
    GPlexQF &outChi2, float *minChi2, int *bestHit, int hit_cnt, int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
#if 0
  getNewBestHitChi2_kernel <<< grid, block, 0, stream >>>
    (outChi2.ptr, minChi2, bestHit, hit_cnt, N);
#endif
}

void fill_array_cu(float *array, int size, int value) {
  thrust::device_ptr<float> d_ptr(array);
  thrust::fill(d_ptr, d_ptr + size, value);
}


__device__ void updateTracksWithBestHit_fn(Hit *hits, int *XHitPos,
    float minChi2, int bestHit,
    float *msErr, int msErr_stride, int msErr_plex_size,
    float *msPar, int msPar_stride, int msPar_plex_size,
    float *propPar, int propPar_stride,
    float *Chi2, int *HitsIdx,
    int N) {
  int itrack = threadIdx.x + blockDim.x*blockIdx.x;
  if (itrack < N) {
    if (bestHit >= 0)
    {
      const char *varr      = (char*) hits;
      const int   off_error = (char*) hits[0].errArrayCU() - varr;
      const int   off_param = (char*) hits[0].posArrayCU() - varr;

      /*Hit   &hit  = hits[ XHitPos.At[itrack] + bestHit[itrack] ];*/
      float &chi2_local = minChi2;
	  
      /*msErr[Nhits].CopyIn(itrack, hit.errArray());*/
      /*SlurpIn_fn<float>(msErr, msErr_stride, msErr_plex_size,
        varr + (itrack*sizeof(Hit)) + off_error, XHitPos, N);*/
      /*msPar[Nhits].CopyIn(itrack, hit.posArray());*/
      /*SlurpIn_fn<float>(msPar, msPar_stride, msPar_plex_size,
        varr + (itrack*sizeof(Hit)) + off_param, XHitPos, N);*/
      for (int i = 0; i < msErr_plex_size; ++i) {
        msErr[i*msErr_stride + itrack] =  hits[XHitPos[itrack]+bestHit].errArrayCU()[i];
      }
      for (int i = 0; i < msPar_plex_size; ++i) {
        msPar[i*msErr_stride + itrack] =  hits[XHitPos[itrack]+bestHit].posArrayCU()[i];
      }
      /*Chi2(itrack, 0, 0) += chi2_local;*/
      Chi2[itrack] += chi2_local;
      /*HitsIdx[Nhits](itrack, 0, 0) = XHitPos.At(itrack, 0, 0) + bestHit[itrack];*/
      HitsIdx[itrack] = XHitPos[itrack] + bestHit;
    }
    else
    {
      /*msErr[Nhits].SetDiagonal3x3(itrack, 666);*/
      msErr[0*msErr_stride + itrack] = 666;
      msErr[1*msErr_stride + itrack] = 0;
      msErr[2*msErr_stride + itrack] = 666;
      msErr[3*msErr_stride + itrack] = 0;
      msErr[4*msErr_stride + itrack] = 0;
      msErr[5*msErr_stride + itrack] = 666;

      /*msPar[Nhits](itrack,0,0) = Par[iP](itrack,0,0);*/
      /*msPar[Nhits](itrack,1,0) = Par[iP](itrack,1,0);*/
      /*msPar[Nhits](itrack,2,0) = Par[iP](itrack,2,0);*/
      for (int i = 0; i < msPar_plex_size; ++i) {
        msPar[i*msPar_stride + itrack] = propPar[i*propPar_stride + itrack]; 
      }
      /*HitsIdx[Nhits](itrack, 0, 0) = -1;*/
      HitsIdx[itrack] = -1;

      // Don't update chi2
    }
  }
}

void updateTracksWithBestHit_wrapper(cudaStream_t &stream,
    BunchOfHitsCU &bunch, GPlexQI &XHitPos, 
    float *minChi2, int *best_hit, 
    GPlexHS &msErr, GPlexHV &msPar,
    GPlexLV &propPar,
    float *Chi2, int *HitsIdx, int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
/*
  updateTracksWithBestHit_kernel <<< grid, block, 0, stream >>>
      (bunch.m_hits, XHitPos.ptr,
       minChi2, best_hit,
       msErr.ptr, msErr.stride, msErr.y,
       msPar.ptr, msPar.stride, msPar.y,
       propPar.ptr, propPar.stride,
       Chi2, HitsIdx,
       N);
*/
}

int getMaxNumHits_wrapper(GPlexQI d_XHitSize, int N) {
  thrust::device_ptr<int> d_ptr(d_XHitSize.ptr);
  int maxSize=  thrust::reduce(d_ptr, d_ptr + N, -1, thrust::maximum<int>());
  maxSize = std::min(maxSize, Config::maxHitsConsidered);

  return maxSize;
}

__global__ void bestHit_kernel(
    Hit *hits, int *XHitPos, 
    float* propErr, size_t propErr_stride,
    float* msErr, size_t msErr_stride, size_t msErr_plex_size,
    float *msPar, size_t msPar_stride, size_t msPar_plex_size,
    float *propPar, size_t propPar_stride,
    float *outChi2, size_t outChi2_stride,
    float *Chi2, int *HitsIdx,
    int maxSize, int N) {

  /*int itrack = threadIdx.x + blockDim.x*blockIdx.x;*/
  int bestHit_reg = -1;
  float minChi2_reg = 15.f;

  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
    HitToMs_fn(msErr, msErr_stride, msErr_plex_size,
               msPar, msPar_stride, msPar_plex_size,
               hits, XHitPos, hit_cnt, N);
#if 0
      // TODO: add CMSGeom
      if (Config::useCMSGeom) {
        //propagateHelixToRMPlex(psErr,  psPar, inChg,  msPar, propErr, propPar);
        throw std::runtime_error("useCMSGeom not implemented yet for GPU");
      } else {}
#endif
    computeChi2_fn(propErr, propErr_stride,
                   msErr, msErr_stride,
                   msPar, msPar_stride,
                   propPar, propPar_stride,
                   outChi2, outChi2_stride,
                   N);
    getNewBestHitChi2_fn(outChi2, minChi2_reg, bestHit_reg, hit_cnt, N);
  }
  updateTracksWithBestHit_fn
      (hits, XHitPos,
       minChi2_reg, bestHit_reg,
       msErr, msErr_stride, msErr_plex_size,
       msPar, msPar_stride, msPar_plex_size,
       propPar, propPar_stride,
       Chi2, HitsIdx,
       N);
}


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
    (bunch.m_hits, XHitPos.ptr,
     propErr.ptr, propErr.stride,
     msErr.ptr, msErr.stride, msErr.kSize,
     msPar.ptr, msPar.stride, msPar.kSize,
     propPar.ptr, propPar.stride,
     outChi2.ptr, outChi2.stride,
     Chi2, HitsIdx,
     maxSize, N);
}

__device__ float downPhi_fn(float phi) {
  while (phi >= Config::PI) {phi-=Config::TwoPI;}
  return phi;
}
	
__device__ float upPhi_fn(float phi) {
  while (phi <= -Config::PI) {phi+=Config::TwoPI;}
  return phi;
}

__device__ float normalizedPhi_fn(float phi) {
  //  return std::fmod(phi, (float) Config::PI); // return phi +pi out of phase for |phi| beyond boundary! 
  if (abs(phi)>=Config::PI) {phi = (phi>0 ? downPhi_fn(phi) : upPhi_fn(phi));}
  return phi;
}

__device__ int getPhiPartition_fn(float phi)
{
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<Config::PI)) std::cout << "anomalous phi=" << phi << std::endl;
  //  const float phiPlusPi  = std::fmod(phi+Config::PI,Config::TwoPI); // normaliztion done here
  const float phiPlusPi = phi+Config::PI; 
  int bin = phiPlusPi*Config::fPhiFactor;
  
  // theoretically these checks below should be taken care of by normalizedPhi, however...
  // these condition checks appeared in very bizarre corner case where propagated phi == pi != Config::PI in check of normalizedPhi (but not unexpected... comparing float point numbers)
  // i.e. delta on floating point check smaller than comparison... making what should be bin = nPhiPart - 1 instead bin = nPhiPart (out of bounds!!) ...or worse if unsigned bin < 0, bin == int max!
  if (bin<0)                      bin = 0;
  else if (bin>=Config::nPhiPart) bin = Config::nPhiPart - 1;

  return bin;
}

__device__ float getPhi_fn(float x, float y)
{
  return atan2(y,x); 
}

__global__ void selectHitRanges_kernel(Hit *hits,
    int *phi_bin_infos_first, int *phi_bin_infos_second, int bunch_fill_index,
    int *XHitPos, int *XHitSize, 
    float *Err, int Err_stride,
    float *Par, int Par_stride,
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

    const float predx = Par[(0*1 + 0)*Par_stride + itrack];  // Par[iI].ConstAt(itrack, 0, 0);
    const float predy = Par[(1*1 + 0)*Par_stride + itrack];  // Par[iI].ConstAt(itrack, 1, 0);
    const float predz = Par[(2*1 + 0)*Par_stride + itrack];  // Par[iI].ConstAt(itrack, 2, 0);

    float phi = getPhi_fn(predx,predy);

    const float px2py2 = predx*predx+predy*predy; // predicted radius^2
    const float dphidx = -predy/px2py2;
    const float dphidy =  predx/px2py2;
    // const float dphi2  =     dphidx*dphidx*(Err[iI].ConstAt(itrack, 0, 0) /*propState.errors.At(0,0)*/) +
    //                          dphidy*dphidy*(Err[iI].ConstAt(itrack, 1, 1) /*propState.errors.At(1,1)*/) +
    //                      2 * dphidx*dphidy*(Err[iI].ConstAt(itrack, 0, 1) /*propState.errors.At(0,1)*/);
    const float dphi2  =     dphidx*dphidx*Err[(0)*Err_stride + itrack] +
                             dphidy*dphidy*Err[(2)*Err_stride + itrack] +
                         2 * dphidx*dphidy*Err[(1)*Err_stride + itrack];

    const float dphi       = sqrtf(fabs(dphi2));//how come I get negative squared errors sometimes? MT -- how small?
    const float nSigmaDphi = fminf(fmaxf(Config::nSigma*dphi, Config::minDPhi), Config::PI);
    //const float nSigmaDphi = Config::nSigma*dphi;

    float dPhiMargin = 0.;
    if (useCMSGeom) {
      //now correct for bending and for layer thickness unsing linear approximation
      /*const float predpx = Par[iP].ConstAt(itrack, 3, 0);*/
      /*const float predpy = Par[iP].ConstAt(itrack, 4, 0);*/
      const float predpx = Par[(3*1 + 0)*Par_stride + itrack];
      const float predpy = Par[(4*1 + 0)*Par_stride + itrack];
      float deltaR = Config::cmsDeltaRad; //fixme! using constant vale, to be taken from layer properties
      float radius = sqrt(px2py2);
      float pt     = sqrt(predpx*predpx + predpy*predpy);
      float cosTheta = ( predx*predpx + predy*predpy )/(pt*radius);
      float hipo = deltaR/cosTheta;
      float dist = sqrt(hipo*hipo - deltaR*deltaR);
      dPhiMargin = dist/radius;
    }
    const float dphiMinus = normalizedPhi_fn(phi-nSigmaDphi-dPhiMargin);
    const float dphiPlus  = normalizedPhi_fn(phi+nSigmaDphi+dPhiMargin);
// FIXME ^ OK

#ifdef DEBUG
    std::ostringstream xout;
    bool               xout_dump = false;
    xout << "--------------------------------------------------------------------------------\n";
    xout << "phi  = " << phi  << ", dphiMinus = " << dphiMinus << ", dphiPlus = " << dphiPlus << std::endl;
    xout << "dphi = " << dphi  << ", dphi2 = " << dphi2 << ", nSigmaDphi = " << nSigmaDphi << ", nSigma = " << Config::nSigma << std::endl;
#endif

    int   phiBinMinus = getPhiPartition_fn(dphiMinus);
    int   phiBinPlus  = getPhiPartition_fn(dphiPlus);

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

void selectHitRanges_wrapper(cudaStream_t &stream, BunchOfHitsCU &bunch, 
    GPlexQI &XHitPos, GPlexQI &XHitSize,
    GPlexLS &Err, GPlexLV &Par,
    int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  selectHitRanges_kernel <<< grid, block, 0, stream >>>
    (bunch.m_hits, bunch.m_phi_bin_infos_first, bunch.m_phi_bin_infos_second, bunch.m_fill_index,
     XHitPos.ptr, XHitSize.ptr, 
     Err.ptr, Err.stride, Par.ptr, Par.stride,
     Config::useCMSGeom, N);
}
