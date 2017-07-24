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

#define BLOCK_SIZE_X 256

using idx_t = Matriplex::idx_t;

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

  int get_Nalloc() const { return Nalloc; }
  void setNumberTracks(const idx_t Ntracks);

  void propagationMerged(const int hit_idx);
  void kalmanUpdateMerged(const int hit_idx);
  void kalmanUpdate_standalone(
      const MPlexLS &psErr, const MPlexLV& psPar, const MPlexQI &inChg,
      const MPlexHS &msErr, const MPlexHV& msPar,
      MPlexLS &outErr, MPlexLV& outPar,
      const int hit_idx, const int N_proc);

  void allocate_extra_addBestHit();
  void free_extra_addBestHit();
  void allocate_extra_combinatorial();
  void free_extra_combinatorial();
  void setHitsIdxToZero(const int hit_idx);

  void addBestHit(EventOfHitsCU& event, GeometryCU &geom_cu,
                  EventOfCandidatesCU &event_of_cands_cu);

  void propagateTracksToR(const float radius, const int N);
  void propagateTracksToR_standalone(const float radius, const int N,
      const MPlexLS& Err_iC, const MPlexLV& par_iC, 
      const MPlexQI& inChg, 
      MPlexLS& Err_iP, MPlexLV& Par_iP);

  // fitting higher order methods
  void FitTracks(Track *tracks_cu, int num_tracks,
                 EventOfHitsCU &events_of_hits_cu,
                 int Nhits, const bool useParamBfield=false);
  void InputTracksAndHitIdx(const EtaBinOfCandidatesCU &etaBin,
                   const int beg, const int end, const bool inputProp);
  void OutputTracksAndHitIdx(EtaBinOfCandidatesCU &etaBin,
                    const int beg, const int end, const bool outputProp);

  void InputTracksAndHitIdxComb(const EtaBinOfCombCandidatesCU &etaBin,
      const int Nhits, const int beg, const int end, const bool inputProp);

  //====================
  // TEMP FCT; copy stuffs
  //====================
  void UpdateWithLastHit(LayerOfHitsCU& layer, int ilay, int N);
  void GetHitFromLayer_standalone(LayerOfHitsCU& layer_cu,
    MPlexQI& HitsIdx, MPlexHV& msPar, MPlexHS& msErr, int hit_idx, int N);
  void UpdateMissingHits_standalone(
      MPlexLS& Err_iP, MPlexLV& Par_iP, 
      MPlexLS& Err_iC, MPlexLV& Par_iC, 
      MPlexQI& HitsIdx, 
      int hit_idx, int N);
  void UpdateWithLastHit_standalone(
      LayerOfHitsCU& layer_cu, MPlexQI& HitsIdx, 
      MPlexLS &Err_iP, MPlexLV& Par_iP, MPlexQI &inChg,
      MPlexHV& msPar, MPlexHS& msErr, 
      MPlexLS &Err_iC, MPlexLV& Par_iC,
      int Nhits, int N_proc);

  void FindTracksInLayers(LayerOfHitsCU *layers, 
                          EventOfCombCandidatesCU& event_of_cands_cu,
                          GeometryCU &geom, bool seed_based=false);
  
  //====================
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
  GPlexHV d_msPar[Config::nLayers];  // on the CPU, with arrays on the GPU
  GPlexHS *d_msErr_arr;
  GPlexHS d_msErr[Config::nLayers];
  
  GPlexQI d_XHitPos;  // QI : 1D arrary following itracks
  GPlexQI d_XHitSize;  // QI : " "
  GPlexHitIdx d_XHitArr;

  GPlexQF d_outChi2;
  GPlexQI *d_HitsIdx_arr;
  GPlexQI d_HitsIdx[Config::nLayers];
  GPlexQF d_Chi2;
  GPlexQI d_Label;

  GPlexQI d_SeedIdx;
  GPlexQI d_CandIdx;

  GPlexQB d_Valid;

  int *d_maxSize;  // max number of tracks for AddBestHit

  // everything run in a stream so multiple instance of FitterCU can
  // run concurrently on the GPU.
  cudaStream_t stream;
};

#include "FitterCU-imp.h"

#endif  // _PROPAGATOR_CU_H_
