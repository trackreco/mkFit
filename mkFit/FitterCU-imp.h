#include <cstdlib>
#include "Config.h"
#include "GeometryCU.h"
#include "reorganize_gplex.h"

template <typename T>
void FitterCU<T>::setNumberTracks(const idx_t Ntracks) {
  N = Ntracks;

  // Raise an exceptioin when the FitterCU instance is too small
  // This should not happen as the loop over tracks in runFittestGPU 
  // takes care of it.
  if (Ntracks > Nalloc) {
    throw std::length_error("FitterCU: Ntracks should be smaller than Nalloc");
  }
}

template <typename T>
void FitterCU<T>::createStream() {
  cudaStreamCreate(&stream);
}

template <typename T>
void FitterCU<T>::destroyStream() {
  cudaStreamDestroy(stream);
}

template <typename T>
void FitterCU<T>::allocateDevice() {
  d_par_iP.allocate(Nalloc, LV);
  d_par_iC.allocate(Nalloc, LV);

  d_Err_iP.allocate(Nalloc, LS);
  d_Err_iC.allocate(Nalloc, LS);

  d_inChg.allocate(Nalloc, QI);
  d_errorProp.allocate(Nalloc, LL);

  cudaMalloc((void**)&d_msPar_arr, Config::nLayers * sizeof(GPlexHV));
  cudaMalloc((void**)&d_msErr_arr, Config::nLayers * sizeof(GPlexHS));
  for (int hi = 0; hi < Config::nLayers; ++hi) {
    d_msPar[hi].allocate(Nalloc, HV);
    d_msErr[hi].allocate(Nalloc, HS);
  }
  cudaMemcpy(d_msPar_arr, d_msPar, Config::nLayers*sizeof(GPlexHV), cudaMemcpyHostToDevice);
  cudaMemcpy(d_msErr_arr, d_msErr, Config::nLayers*sizeof(GPlexHS), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_maxSize, sizeof(int));  //  global maximum
  cudaCheckError();
}

template <typename T>
void FitterCU<T>::freeDevice() {
  d_par_iC.free();
  d_inChg.free();
  d_par_iP.free();
  d_errorProp.free();
  d_Err_iP.free();
  d_Err_iC.free();
  for (int hi = 0; hi < Config::nLayers; ++hi) {
    d_msPar[hi].free();
    d_msErr[hi].free();
  }
  cudaFree(d_msPar_arr);
  cudaFree(d_msErr_arr);
  cudaFree(d_maxSize);
  cudaCheckError();
}

template <typename T>
void FitterCU<T>::kalmanUpdateMerged(const int hit_idx) {
  kalmanUpdate_wrapper(stream, d_Err_iP, d_msErr[hit_idx],
                       d_par_iP, d_msPar[hit_idx], d_par_iC, d_Err_iC, N);
}

template <typename T>
void FitterCU<T>::kalmanUpdate_standalone(
    const MPlexLS &psErr, const MPlexLV& psPar, const MPlexQI &inChg,
    const MPlexHS &msErr, const MPlexHV& msPar,
    MPlexLS &outErr, MPlexLV& outPar,
    const int hit_idx, const int N_proc)
{
  d_Err_iP.copyAsyncFromHost(stream, psErr);
  d_msErr[hit_idx].copyAsyncFromHost(stream, msErr);
  d_par_iP.copyAsyncFromHost(stream, psPar);
  d_msPar[hit_idx].copyAsyncFromHost(stream, msPar);

  kalmanUpdate_wrapper(stream, d_Err_iP, d_msErr[hit_idx],
                       d_par_iP, d_msPar[hit_idx], d_par_iC, d_Err_iC, N_proc);

  d_par_iC.copyAsyncToHost(stream, outPar);
  d_Err_iC.copyAsyncToHost(stream, outErr);
}

template <typename T>
void FitterCU<T>::propagationMerged(const int hit_idx) {
  propagation_wrapper(stream, d_msPar[hit_idx], d_par_iC, d_inChg,
                      //d_par_iP, d_Err_iC, d_Err_iP, N); // TODO: Check outErr/errorProp
                      d_par_iP, d_errorProp, d_Err_iP, N);
}

// FIXME: Temporary. Separate allocations / transfers
template <typename T>
void FitterCU<T>::allocate_extra_addBestHit() {
  d_outChi2.allocate(Nalloc, QF);
  d_XHitPos.allocate(Nalloc, QI);
  d_XHitSize.allocate(Nalloc, QI);
  d_XHitArr.allocate(Nalloc, GPlexHitIdxMax);
  cudaMalloc((void**)&d_HitsIdx_arr, Config::nLayers * sizeof(GPlexQI));
  for (int hi = 0; hi < Config::nLayers; ++hi) {
    d_HitsIdx[hi].allocate(Nalloc, QI);
  }
  cudaMemcpy(d_HitsIdx_arr, d_HitsIdx, Config::nLayers*sizeof(GPlexQI), cudaMemcpyHostToDevice);
  d_Chi2.allocate(Nalloc, QF);
  d_Label.allocate(Nalloc, QI);
  cudaCheckError();
}

template <typename T>
void FitterCU<T>::free_extra_addBestHit() {
  for (int hi = 0; hi < Config::nLayers; ++hi) {
    d_HitsIdx[hi].free(); cudaCheckError();
  }
  cudaFree(d_HitsIdx_arr);
  d_Label.free(); cudaCheckError();
  d_Chi2.free(); cudaCheckError();

  d_XHitArr.free(); cudaCheckError();
  d_XHitSize.free(); cudaCheckError();
  d_XHitPos.free(); cudaCheckError();
  d_outChi2.free(); cudaCheckError();
}

template <typename T>
void FitterCU<T>::setHitsIdxToZero(const int hit_idx) {
  cudaMemset(d_HitsIdx[hit_idx].ptr, 0, Nalloc*sizeof(int));
}

template <typename T>
void FitterCU<T>::addBestHit(EventOfHitsCU &event, GeometryCU &geom_cu,
                             EventOfCandidatesCU &event_of_cands_cu) {
    findBestHit_wrapper(stream, event.m_layers_of_hits, 
                        event_of_cands_cu,
                        d_XHitSize, d_XHitArr,
                        d_Err_iP, d_par_iP, 
                        d_msErr_arr, d_msPar_arr,
                        d_Err_iC, d_par_iC, d_outChi2,
                        d_Chi2, d_HitsIdx_arr,
                        d_inChg, d_Label, geom_cu,
                        d_maxSize, N);
}   

template <typename T>
void FitterCU<T>::InputTracksAndHitIdx(const EtaBinOfCandidatesCU &etaBin,
                              const int beg, const int end, const bool inputProp) {
  InputTracksCU_wrapper(stream, etaBin, d_Err_iP, d_par_iP,
                        d_inChg, d_Chi2, d_Label, d_HitsIdx_arr,
                        beg, end, inputProp, N);
}

template <typename T>
void FitterCU<T>::OutputTracksAndHitIdx(EtaBinOfCandidatesCU &etaBin,
                               const int beg, const int end, const bool outputProp) {
  OutputTracksCU_wrapper(stream, etaBin, d_Err_iP, d_par_iP,
                         d_inChg, d_Chi2, d_Label, d_HitsIdx_arr,
                         beg, end, outputProp, N);
  cudaStreamSynchronize(stream);
  cudaCheckError();
}


#if 1
template <typename T>
void FitterCU<T>::propagateTracksToR(const float radius, const int N) {
  propagationForBuilding_wrapper(stream, d_Err_iC, d_par_iC, d_inChg, 
                                 radius, d_Err_iP, d_par_iP, N); 
}
#endif

template <typename T>
void FitterCU<T>::propagateTracksToR_standalone(const float radius, const int N,
    const MPlexLS& Err_iC, const MPlexLV& par_iC, const MPlexQI& inChg, 
    MPlexLS& Err_iP, MPlexLV& Par_iP) {
  d_Err_iC.copyAsyncFromHost(stream, Err_iC);
  d_par_iC.copyAsyncFromHost(stream, par_iC);
  //propagationForBuilding_wrapper(stream, d_Err_iC, d_par_iC, d_inChg, 
                                 //radius, d_Err_iP, d_par_iP, N); 
  d_Err_iP.copyAsyncToHost(stream, Err_iP);
  d_par_iP.copyAsyncToHost(stream, Par_iP);
}

template <typename T>
void FitterCU<T>::FitTracks(MPlexQI &Chg, MPlexLV& par_iC, MPlexLS& err_iC,
                            MPlexHV* msPar, MPlexHS* msErr, int Nhits,
                            std::vector<Track> &tracks, int beg, int end,
                            std::vector<HitVec> &layerHits) {
  float etime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  //launch everything in a stream to enable concurrent execution of events
  createStream();

  // allocateDevice();  -> moved to mkFit/mkFit.cc

  setNumberTracks(end-beg);

  d_inChg.copyAsyncFromHost(stream, Chg);
  d_par_iC.copyAsyncFromHost(stream, par_iC);
  d_Err_iC.copyAsyncFromHost(stream, err_iC);

  cudaEventRecord(start, 0);
 
  double total_reorg = 0.;
  for (int hi = 0; hi < Nhits; ++hi)
  {
    // Switch outPut and inPut parameters and errors
    // similar to iC <-> iP in the CPU code.
    d_par_iP.copyAsyncFromDevice(stream, d_par_iC); 
    d_Err_iP.copyAsyncFromDevice(stream, d_Err_iC);
    
    double time_input = dtime();
    int itrack;
    omp_set_num_threads(Config::numThreadsReorg);
#pragma omp parallel for
    for (int i = beg; i < end; ++i) {
      itrack = i - beg;
      Track &trk = tracks[i];

      const int hidx = trk.getHitIdx(hi);
      const Hit &hit = layerHits[hi][hidx];

      msErr[hi].CopyIn(itrack, hit.errArray());
      msPar[hi].CopyIn(itrack, hit.posArray());
    }
    total_reorg += (dtime() - time_input)*1e3;

    d_msPar[hi].copyAsyncFromHost(stream, msPar[hi]);
    d_msErr[hi].copyAsyncFromHost(stream, msErr[hi]);

    propagationMerged(hi);
    kalmanUpdateMerged(hi);
  }
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  cudaEventElapsedTime(&etime, start, stop);
  std::cerr << "CUDA etime: " << etime << " ms.\n";
  std::cerr << "Total reorg: " << total_reorg << " ms.\n";

  d_par_iC.copyAsyncToHost(stream, par_iC);
  d_Err_iC.copyAsyncToHost(stream, err_iC);
  
  cudaStreamSynchronize(stream);
  // freeDevice(); -> moved to mkFit/mkFit.cc
  destroyStream();

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}

///////////////////////////////////////////////////////////////////////////////
// Backup function: temps that have been deactivated
///////////////////////////////////////////////////////////////////////////////

#if 0
template <typename T>
void FitterCU<T>::computeChi2gpu(const MPlexLS &psErr, const MPlexLV& propPar,
    const MPlexQI &inChg, MPlexHS &msErr, MPlexHV& msPar,
    float *minChi2, int *bestHit,
    LayerOfHitsCU &d_layer, MPlexQI &XHitSize, Matriplex::Matriplex<int, 16, 1, MPT_SIZE> &XHitArr,
    MPlexQF &Chi2, MPlexQI &HitsIdx, MPlexQF &outChi2, int maxSize2, int hit_idx, int NN) {

  float *d_minChi2;
  int *d_bestHit;
  cudaMalloc((void**)&d_minChi2, NN*sizeof(float));
  cudaMalloc((void**)&d_bestHit, NN*sizeof(int));

  cudaMemcpyAsync(d_minChi2, minChi2, NN*sizeof(float), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(d_bestHit, bestHit, NN*sizeof(int), cudaMemcpyHostToDevice, stream);

  cudaMemset(d_bestHit, -1, NN*sizeof(int));
  fill_array_cu(d_minChi2, NN, 15.f);

  d_Err_iP.copyAsyncFromHost(stream, psErr);
  d_par_iP.copyAsyncFromHost(stream, propPar);
  d_msErr[hit_idx].copyAsyncFromHost(stream, msErr);
  d_msPar[hit_idx].copyAsyncFromHost(stream, msPar);
  //d_XHitPos.copyAsyncFromHost(stream, XHitPos);
  //d_XHitSize.copyAsyncFromHost(stream, XHitSize);
  //d_XHitArr.copyAsyncFromHost(stream, XHitArr);

  //cudaMemcpy2DAsync(d_Chi2, NN*sizeof(float), Chi2.fArray, NN*sizeof(float), 
               //NN*sizeof(float), 1, cudaMemcpyHostToDevice, stream);
  //cudaMemcpy2DAsync(d_HitsIdx, NN*sizeof(int), HitsIdx.fArray, NN*sizeof(int), 
               //NN*sizeof(int), 1, cudaMemcpyHostToDevice, stream);

  //cudaStreamSynchronize(stream);
  //cudaCheckError();

  //selectHitRanges_wrapper(stream, d_bunch, d_XHitPos, d_XHitSize, 
      //d_Err_iP, d_par_iP, N);

  int maxSize = getMaxNumHits_wrapper(d_XHitSize, N);
  //bestHit_wrapper(stream, d_bunch, d_XHitPos,
                  //d_Err_iP, d_msErr, d_msPar, d_par_iP, d_outChi2,
                  //d_Chi2, d_HitsIdx,
                  //maxSize2, N);
  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
    //// TODO: add CMSGeom
    //if (Config::useCMSGeom) {
      ////propagateHelixToRMPlex(psErr,  psPar, inChg,  msPar, propErr, propPar);
      //throw std::runtime_error("useCMSGeom not implemented yet for GPU");
    //} else {}
    HitToMs_wrapper(stream, d_msErr[hit_idx], d_msPar[hit_idx], d_layer,
                    d_XHitSize, d_XHitArr, d_HitsIdx[hit_idx], hit_cnt, NN);

    computeChi2_wrapper(stream, d_Err_iP, d_msErr[hit_idx], //d_resErr, 
        d_msPar[hit_idx], d_par_iP, d_outChi2, NN);

    getNewBestHitChi2_wrapper(stream, d_XHitSize, d_XHitArr, d_outChi2, d_minChi2, d_bestHit, hit_cnt, NN);

    //cudaStreamSynchronize(stream);
    //cudaCheckError();
  }
  updateTracksWithBestHit_wrapper(stream, d_layer, d_minChi2, d_bestHit, 
    d_msErr[hit_idx], d_msPar[hit_idx], d_par_iP, d_Chi2, d_HitsIdx[hit_idx], N);

  //d_outChi2.copyAsyncToHost(stream, outChi2);
  //cudaMemcpyAsync(minChi2, d_minChi2, NN*sizeof(float), cudaMemcpyDeviceToHost, stream);
  //cudaMemcpyAsync(bestHit, d_bestHit, NN*sizeof(int), cudaMemcpyDeviceToHost, stream);

  //cudaMemcpy2DAsync(Chi2.fArray, NN*sizeof(float), d_Chi2, NN*sizeof(float), 
  //             NN*sizeof(float), 1, cudaMemcpyDeviceToHost, stream);
  //cudaMemcpy2DAsync(HitsIdx.fArray, NN*sizeof(int), d_HitsIdx, NN*sizeof(int), 
  //             NN*sizeof(int), 1, cudaMemcpyDeviceToHost, stream);
  d_Chi2.copyAsyncToHost(stream, Chi2); 
  d_HitsIdx[hit_idx].copyAsyncToHost(stream, HitsIdx);
  d_msErr[hit_idx].copyAsyncToHost(stream, msErr);
  d_msPar[hit_idx].copyAsyncToHost(stream, msPar);


  cudaStreamSynchronize(stream);
  cudaCheckError();
  //for (int itrack = 0; itrack < NN; ++itrack)
  //{
    ////printf("CPU [%d]  -- %d : %f\n", itrack, HitsIdx(itrack, 0, 0), Chi2[itrack]);
  //}

  cudaFree(d_minChi2);
  cudaFree(d_bestHit);
}

// FIXME: Temporary. Separate allocations / transfers
template <typename T>
void FitterCU<T>::prepare_addBestHit() {
    //const MPlexLS &psErr, const MPlexLV& propPar,
    //const MPlexQI &inChg,
    //MPlexQI &XHitSize, Matriplex::Matriplex<int, 16, 1, MPT_SIZE> &XHitArr,
    //size_t num_tracks) {
  //setNumberTracks(num_tracks);  // temporary: should be end - beg

  //createStream();
  //cudaCheckError();
  // psErr -> d_Err_iP
  //d_Err_iP.copyAsyncFromHost(stream, psErr);
  //d_par_iP.copyAsyncFromHost(stream, propPar);
  //d_inChg.copyAsyncFromHost(stream, inChg);
}

// TODO: Temporary. Separate allocations / transfers
template <typename T>
void FitterCU<T>::finalize_addBestHit(
    MPlexHS *msErr, MPlexHV* msPar,
    MPlexLS& Err_iP, MPlexLV& Par_iP, 
    MPlexQI *HitsIdx, 
    MPlexQI &Label,
    int start_idx, int end_idx) {
  d_par_iP.copyAsyncToHost(stream, Par_iP);
  d_Err_iP.copyAsyncToHost(stream, Err_iP);
  d_Label.copyAsyncToHost(stream, Label);
 
  // Get msPar, msErr, chi2 and HitIdx out from the GPU to the CPU
  for (int hit_idx = start_idx; hit_idx < end_idx; ++hit_idx) {
    d_msPar[hit_idx].copyAsyncToHost(stream, msPar[hit_idx]);
    d_msErr[hit_idx].copyAsyncToHost(stream, msErr[hit_idx]);
    d_HitsIdx[hit_idx].copyAsyncToHost(stream, HitsIdx[hit_idx]);
  }
}
#endif
