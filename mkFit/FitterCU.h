#ifndef _PROPAGATOR_CU_H_
#define _PROPAGATOR_CU_H_

#include <cuda_runtime.h>

#include "Matrix.h"

#include "HitStructuresCU.h"
#include "GPlex.h"
#include "GeometryCU.h"
#include "gpu_utils.h"

#include "propagation_kernels.h"
#include "kalmanUpdater_kernels.h"
#include "computeChi2_kernels.h"
#include "index_selection_kernels.h"
#include "best_hit_kernels.h"

#include <omp.h>
#include <stdexcept>

constexpr int LV = 6;
constexpr int QI = 1;
constexpr int QF = 1;
#define LL 36
constexpr int LS = 21;
constexpr int HV = 3;
constexpr int HS = 6;
constexpr int LH = 18;

#define BLOCK_SIZE_X 256

using idx_t = Matriplex::idx_t;

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
  cudaStream_t& get_stream() { return stream; }

  void setNumberTracks(const idx_t Ntracks);

  void propagationMerged(const int hit_idx);
  void kalmanUpdateMerged(const int hit_idx);
  void kalmanUpdate_standalone(
      const MPlexLS &psErr, const MPlexLV& psPar, const MPlexQI &inChg,
      const MPlexHS &msErr, const MPlexHV& msPar,
      MPlexLS &outErr, MPlexLV& outPar,
      const int hit_idx, const int N_proc);

#if 0
  void computeChi2gpu(const MPlexLS &psErr, const MPlexLV& propPar,
    const MPlexQI &inChg, MPlexHS &msErr, MPlexHV& msPar,
    float *minChi2, int *bestHit,
    LayerOfHitsCU &d_layer, MPlexQI &XHitSize, Matriplex::Matriplex<int, 16, 1, MPT_SIZE> &XHitArr,
    MPlexQF &Chi2, MPlexQI &HitsIdx, MPlexQF&outChi2, int maxSize, int hit_idx,
    int NN);
#endif

  void allocate_extra_addBestHit();
  void free_extra_addBestHit();

#if 0
  void prepare_addBestHit();
      //const MPlexLS &psErr, const MPlexLV& propPar,
      //const MPlexQI &inChg, 
      //MPlexQI &XHitSize, Matriplex::Matriplex<int, 16, 1, MPT_SIZE> &XHitArr,
      //size_t NN);
  void finalize_addBestHit(
      MPlexHS *msErr, MPlexHV* msPar,
      MPlexLS& Err_iP, MPlexLV& Par_iP, 
      MPlexQI *HitsIdx,
      MPlexQI &Label,
      int start_idx, int end_idx);
#endif
  void setHitsIdxToZero(const int hit_idx);

#if 1
  //void addBestHit(EventOfHitsCU& event, const int ilay, const float *radii, int hit_idx);
  void addBestHit(EventOfHitsCU& event, GeometryCU &geom_cu,
                  EventOfCandidatesCU &event_of_cands_cu);
#endif
  void propagateTracksToR(const float radius, const int N);
  void propagateTracksToR_standalone(const float radius, const int N,
      const MPlexLS& Err_iC, const MPlexLV& par_iC, 
      const MPlexQI& inChg, 
      MPlexLS& Err_iP, MPlexLV& Par_iP);

  // fitting higher order methods
  void FitTracks(MPlexQI &Chg, MPlexLV& par_iC, MPlexLS& err_iC,
                 MPlexHV* msPar, MPlexHS* msErr, int Nhits,
                 std::vector<Track> &tracks, int beg, int end,
                 std::vector<HitVec> &layerHits);
  void InputTracksAndHitIdx(const EtaBinOfCandidatesCU &etaBin,
                   const int beg, const int end, const bool inputProp);
  void OutputTracksAndHitIdx(EtaBinOfCandidatesCU &etaBin,
                    const int beg, const int end, const bool outputProp);

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

  GPlexHV *d_msPar_arr;  // completely on the GPU
  GPlexHV d_msPar[MAX_HITS];  // on the CPU, with arrays on the GPU
  GPlexHS *d_msErr_arr;
  GPlexHS d_msErr[MAX_HITS];
  
  GPlexQI d_XHitPos;  // QI : 1D arrary following itracks
  GPlexQI d_XHitSize;  // QI : " "
  GPlexHitIdx d_XHitArr;

  GPlexQF d_outChi2;
  GPlexQI *d_HitsIdx_arr;
  GPlexQI d_HitsIdx[MAX_HITS];
  GPlexQF d_Chi2;
  GPlexQI d_Label;

  int *d_maxSize;  // max number of tracks for AddBestHit

  // everything run in a stream so multiple instance of FitterCU can
  // run concurrently on the GPU.
  cudaStream_t stream;
};

#include "FitterCU-imp.h"

#endif  // _PROPAGATOR_CU_H_
