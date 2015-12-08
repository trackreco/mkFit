#ifndef _PROPAGATOR_CU_H_
#define _PROPAGATOR_CU_H_

#include <cuda_runtime.h>
#include <omp.h>

#include "Matrix.h"
#include "propagation_kernels.h"
#include "kalmanUpdater_kernels.h"
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

//Macro for checking cuda errors following a cuda launch or api call
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

// FIXME: remove me
//void raiseNotImplemented() {
  //throw ("notImplemented error"); 
//}


template <typename T>
class FitterCU {
 public:
  FitterCU (idx_t N) : N(N) {};
  virtual ~FitterCU () {};

  void allocateDevice();
  void freeDevice();

  void createStream();
  void destroyStream();

  void sendInParToDevice(const MPlexLV& inPar);
  void sendInErrToDevice(const MPlexLS& inErr);
  void sendInChgToDevice(const MPlexQI& inChg);
  void sendMsRadToDevice(const MPlexQF& msRad);
  void sendOutParToDevice(const MPlexLV& outPar); // TODO: temporary
  // TODO: actually inErr should be send to the device as initial
  //       values for outErr. Rename accordingly?
  void sendOutErrToDevice(const MPlexLS& outErr);
  void sendMsParToDevice(const MPlexHV& msPar);
  
  //void getOutParFromDevice(MPlexLV& outPar);  // Par_iP 
  void getErrorPropFromDevice(MPlexLL& errorProp);
  //void getOutErrFromDevice(MPlexLS& outErr);  // Err_iP
  void getMsRadFromDevice(MPlexQF& msRad);

  void setOutParFromInPar();
  void setOutErrFromInErr();

  // updater specfic transfers.
  void sendMsErrToDevice(const MPlexHS& msErr);
  void getOutParFromDevice(MPlexLV& outPar);
  void getOutErrFromDevice(MPlexLS& outErr);

  void swap_iP_iC();

  // Propagator methods
  void computeMsRad();
  void helixAtRFromIterative();
  void similarity();

  // Updater methods
#if 0
  void addIntoUpperLeft3x3();
  void InvertCramerSym();
  void multKalmanGainCU();
  void multResidualsAdd();
  void kalmanGain_x_propErr();
#endif

  void propagationMerged();
  void kalmanUpdateMerged();

  // fitting higher order methods
  void FitTracks(MPlexQI &Chg, MPlexLV& par_iC, MPlexLS& err_iC,
                 MPlexHV* msPar, MPlexHS* msErr, int Nhits,
                 std::vector<Track> &tracks, int beg, int end,
                 std::vector<HitVec> &layerHits);

 private:
  idx_t N;
  /* data */
  GPlex<T> d_par_iC;  // LV
  GPlex<int> d_inChg;  // QI
  GPlex<T> d_par_iP; // LV
  GPlex<T> d_msRad;  // QF
  GPlex<T> d_errorProp;  // LL
  GPlex<T> d_Err_iP;
  GPlex<T> d_msPar;

  // variables from KalmanUpdaterCU.h not shared with FitterCU
  GPlex<T> d_kalmanGain;
  GPlex<T> d_outErr;
  GPlex<T> d_resErr;
  GPlex<T> d_msErr;
  
  // everything run in a stream so multiple instance of FitterCU can
  // run concurrently on the GPU.
  cudaStream_t stream;
};

#include "FitterCU-imp.h"

#endif  // _PROPAGATOR_CU_H_
