#ifndef _PROPAGATOR_CU_H_
#define _PROPAGATOR_CU_H_

#include <cuda_runtime.h>
#include <omp.h>
#include <stdexcept>

#include "Matrix.h"
#include "propagation_kernels.h"
#include "kalmanUpdater_kernels.h"
#include "computeChi2_kernels.h"
#include "HitStructuresCU.h"
#include "GPlex.h"

#define LV 6
#define QI 1
#define QF 1
#define LL 36
#define LS 21
#define HV 3
#define HS 6
#define LH 18

#define BLOCK_SIZE_X 16
#define MAX_BLOCKS_X 65535 // CUDA constraint

using idx_t = Matriplex::idx_t;

// Macro for checking cuda errors following a cuda launch or api call
// This comes from Jeff Larkins (NVIDIA)
#define cudaCheckError() {                                          \
  cudaError_t e=cudaGetLastError();                                 \
  if(e!=cudaSuccess) {                                              \
    printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
    exit(0); \
  }                                                                 \
}
#if 0
#define cudaCheckErrorSync() {                                      \
  cudaDeviceSynchronize();                                          \
  cudaCheckError();                                                 \
}
#else
#define cudaCheckErrorSync() {}
#endif

void separate_first_call_for_meaningful_profiling_numbers();

template <typename T>
class FitterCU {
 public:
  FitterCU (idx_t Nalloc) : Nalloc(Nalloc) {};
  virtual ~FitterCU () {};

  void allocateDevice();
  void freeDevice();

  void createStream();
  void destroyStream();

  void setNumberTracks(idx_t Ntracks);

  void propagationMerged();
  void kalmanUpdateMerged();
  void kalmanUpdate_standalone(
      const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
      const MPlexHS &msErr,  const MPlexHV& msPar,
      MPlexLS &outErr,       MPlexLV& outPar);

  void computeChi2gpu(const MPlexLS &psErr, const MPlexLV& propPar,
    const MPlexQI &inChg, MPlexHS &msErr, MPlexHV& msPar,
    BunchOfHitsCU &d_bunch, //MPlexQI &XHitPos, MPlexQI &XHitSize,
    MPlexQF &Chi2, MPlexQI &HitsIdx,
    int NN);

  void allocate_extra_addBestHit();
  void free_extra_addBestHit();

  void prepare_addBestHit(
      const MPlexLS &psErr, const MPlexLV& propPar,
      const MPlexQI &inChg, 
      size_t NN);
  void finalize_addBestHit(
      MPlexHS &msErr, MPlexHV& msPar,
      MPlexLS &outErr, MPlexLV &outPar,
      MPlexQI &HitsIdx, MPlexQF &Chi2);

  void addBestHit(BunchOfHitsCU &bunch);

  // fitting higher order methods
  void FitTracks(MPlexQI &Chg, MPlexLV& par_iC, MPlexLS& err_iC,
                 MPlexHV* msPar, MPlexHS* msErr, int Nhits,
                 std::vector<Track> &tracks, int beg, int end,
                 std::vector<HitVec> &layerHits);

 private:
  // N is the actual size, Nalloc should be >= N, as it is intended
  // to allocated arrays that can be used for several sets of tracks.
  idx_t Nalloc;
  idx_t N;

  /* data */
  GPlexLV d_par_iP; // LV
  GPlexLV d_par_iC; // LV

  GPlexLS d_Err_iP; // LS
  GPlexLS d_Err_iC; // LS

  GPlexQI d_inChg;  // QI
  GPlexQF d_msRad;  // QF
  GPlexLL d_errorProp;  // LL

  GPlexHV d_msPar;
  GPlexHS d_msErr;
  
  GPlexQI d_XHitPos;  // QI : 1D arrary following itracks
  GPlexQI d_XHitSize;  // QI : " "
  GPlexQF d_outChi2;
  int *d_HitsIdx;
  float *d_Chi2;

  // everything run in a stream so multiple instance of FitterCU can
  // run concurrently on the GPU.
  cudaStream_t stream;
};

#include "FitterCU-imp.h"

#endif  // _PROPAGATOR_CU_H_
