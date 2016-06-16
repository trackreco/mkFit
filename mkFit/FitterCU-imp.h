#include <cstdlib>
#include "Config.h"

template <typename T>
void FitterCU<T>::setNumberTracks(idx_t Ntracks) {
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
  d_msPar.allocate(Nalloc, HV);
  d_msErr.allocate(Nalloc, HS);

  cudaCheckError()
}

template <typename T>
void FitterCU<T>::freeDevice() {
  d_par_iC.free();
  d_inChg.free();
  d_par_iP.free();
  d_errorProp.free();
  d_Err_iP.free();
  d_msPar.free();
  d_Err_iC.free();
  d_msErr.free();

  cudaCheckError()
}

template <typename T>
void FitterCU<T>::kalmanUpdateMerged() {
  kalmanUpdate_wrapper(stream, d_Err_iP, d_msErr,
                       d_par_iP, d_msPar, d_par_iC, d_Err_iC, N);
}

template <typename T>
void FitterCU<T>::kalmanUpdate_standalone(
    const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
    const MPlexHS &msErr,  const MPlexHV& msPar,
    MPlexLS &outErr,       MPlexLV& outPar, int N_proc)
{
  d_Err_iP.copyAsyncFromHost(stream, psErr);
  d_msErr.copyAsyncFromHost(stream, msErr);
  d_par_iP.copyAsyncFromHost(stream, psPar);
  d_msPar.copyAsyncFromHost(stream, msPar);

  kalmanUpdate_wrapper(stream, d_Err_iP, d_msErr,
                       d_par_iP, d_msPar, d_par_iC, d_Err_iC, N_proc);

  d_par_iC.copyAsyncToHost(stream, outPar);
  d_Err_iC.copyAsyncToHost(stream, outErr);
}

template <typename T>
void FitterCU<T>::propagationMerged() {
  propagation_wrapper(stream, d_msPar, d_par_iC, d_inChg,
                      //d_par_iP, d_Err_iC, d_Err_iP, N); // TODO: Check outErr/errorProp
                      d_par_iP, d_errorProp, d_Err_iP, N);
}

#if 1
template <typename T>
void FitterCU<T>::computeChi2gpu(const MPlexLS &psErr, const MPlexLV& propPar,
    const MPlexQI &inChg, MPlexHS &msErr, MPlexHV& msPar,
    float *minChi2, int *bestHit,
    LayerOfHitsCU &d_layer, MPlexQI &XHitSize, Matriplex::Matriplex<int, 16, 1, MPT_SIZE> &XHitArr,
    MPlexQF &Chi2, MPlexQI &HitsIdx, MPlexQF &outChi2, int maxSize2, int NN) {

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
  d_msErr.copyAsyncFromHost(stream, msErr);
  d_msPar.copyAsyncFromHost(stream, msPar);
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
    HitToMs_wrapper(stream, d_msErr, d_msPar, d_layer, d_XHitSize, d_XHitArr, d_HitsIdx, hit_cnt, NN);

    computeChi2_wrapper(stream, d_Err_iP, d_msErr, //d_resErr, 
        d_msPar, d_par_iP, d_outChi2, NN);

    getNewBestHitChi2_wrapper(stream, d_XHitSize, d_XHitArr, d_outChi2, d_minChi2, d_bestHit, hit_cnt, NN);

    //cudaStreamSynchronize(stream);
    //cudaCheckError();
  }
  updateTracksWithBestHit_wrapper(stream, d_layer, d_minChi2, d_bestHit, 
    d_msErr, d_msPar, d_par_iP, d_Chi2, d_HitsIdx, N);

  //d_outChi2.copyAsyncToHost(stream, outChi2);
  //cudaMemcpyAsync(minChi2, d_minChi2, NN*sizeof(float), cudaMemcpyDeviceToHost, stream);
  //cudaMemcpyAsync(bestHit, d_bestHit, NN*sizeof(int), cudaMemcpyDeviceToHost, stream);

  cudaMemcpy2DAsync(Chi2.fArray, NN*sizeof(float), d_Chi2, NN*sizeof(float), 
               NN*sizeof(float), 1, cudaMemcpyDeviceToHost, stream);
  cudaMemcpy2DAsync(HitsIdx.fArray, NN*sizeof(int), d_HitsIdx, NN*sizeof(int), 
               NN*sizeof(int), 1, cudaMemcpyDeviceToHost, stream);
  d_msErr.copyAsyncToHost(stream, msErr);
  d_msPar.copyAsyncToHost(stream, msPar);


  cudaStreamSynchronize(stream);
  cudaCheckError();
  //for (int itrack = 0; itrack < NN; ++itrack)
  //{
    ////printf("CPU [%d]  -- %d : %f\n", itrack, HitsIdx(itrack, 0, 0), Chi2[itrack]);
  //}

  cudaFree(d_minChi2);
  cudaFree(d_bestHit);
}
#endif

// FIXME: Temporary. Separate allocations / transfers
template <typename T>
void FitterCU<T>::allocate_extra_addBestHit() {
  d_outChi2.allocate(Nalloc, QF);
  d_XHitPos.allocate(Nalloc, QI);
  d_XHitSize.allocate(Nalloc, QI);
  d_XHitArr.allocate(Nalloc, GPlexHitIdxMax);
  // FIXME: Make those GPlex-es. and use .allocate()
  cudaMalloc((void**)&d_HitsIdx, Nalloc*sizeof(int)); cudaCheckError();
  cudaMalloc((void**)&d_Chi2, Nalloc*sizeof(float)); cudaCheckError();
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::free_extra_addBestHit() {
  destroyStream();

  cudaFree(d_HitsIdx); cudaCheckError();
  cudaFree(d_Chi2); cudaCheckError();

  d_XHitArr.free(); cudaCheckError();
  d_XHitSize.free(); cudaCheckError();
  d_XHitPos.free(); cudaCheckError();
  d_outChi2.free(); cudaCheckError();
}

// FIXME: Temporary. Separate allocations / transfers
template <typename T>
void FitterCU<T>::prepare_addBestHit(
    const MPlexLS &psErr, const MPlexLV& propPar,
    const MPlexQI &inChg,
    MPlexQI &XHitSize, Matriplex::Matriplex<int, 16, 1, MPT_SIZE> &XHitArr,
    size_t NN) {
  setNumberTracks(NN);  // temporary: should be end - beg

  createStream();
  cudaCheckError()
#if 1
  // psErr -> d_Err_iP
  //cudaMemcpy2DAsync(d_Err_iP.ptr, d_Err_iP.pitch, psErr.fArray, N*sizeof(T),
               //N*sizeof(T), LS, cudaMemcpyHostToDevice, stream);
  d_Err_iP.copyAsyncFromHost(stream, psErr);
  d_par_iP.copyAsyncFromHost(stream, propPar);
  d_inChg.copyAsyncFromHost(stream, inChg);
  
  //cudaMemset2D(d_XHitSize.ptr, d_XHitSize.pitch,
               //0, sizeof(int)*d_XHitSize.N, d_XHitSize.kSize);
  //d_XHitSize.copyAsyncFromHost(stream, XHitSize);
  //d_XHitArr.copyAsyncFromHost(stream, XHitArr);
#endif
}

// TODO: Temporary. Separate allocations / transfers
template <typename T>
void FitterCU<T>::finalize_addBestHit(
    MPlexHS &msErr, MPlexHV& msPar,
    MPlexLS& Err_iC, MPlexLV& Par_iC, 
    MPlexLS& Err_iP, MPlexLV& Par_iP, 
    MPlexQI &HitsIdx, MPlexQF &Chi2) {
#if 1
  d_par_iC.copyAsyncToHost(stream, Par_iC);
  d_Err_iC.copyAsyncToHost(stream, Err_iC);

  d_par_iP.copyAsyncToHost(stream, Par_iP);
  d_Err_iP.copyAsyncToHost(stream, Err_iP);
 
  // Get msPar, msErr, chi2 and HitIdx out from the GPU to the CPU
  d_msPar.copyAsyncToHost(stream, msPar);
  d_msErr.copyAsyncToHost(stream, msErr);
  cudaMemcpyAsync(HitsIdx.fArray, d_HitsIdx, N*sizeof(int), cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(Chi2.fArray, d_Chi2, N*sizeof(float), cudaMemcpyDeviceToHost, stream);
#endif
}

template <typename T>
void FitterCU<T>::setHitsIdxToZero() {
  cudaMemset(d_HitsIdx, 0, Nalloc*sizeof(int));
}

#if 1
template <typename T>
void FitterCU<T>::addBestHit(LayerOfHitsCU &layer, const int ilay, const float radius) {

  //selectHitRanges_wrapper(stream, layer, d_XHitPos, d_XHitSize, 
      //d_Err_iP, d_par_iP, N);
  selectHitIndices_wrapper(stream, 
      layer, d_Err_iP, d_par_iP, 
      d_XHitSize, d_XHitArr, N);

  // TODO: get this thing inside bestHit_kernel
  int maxSize = getMaxNumHits_wrapper(d_XHitSize, N);
  cudaDeviceSynchronize(); cudaCheckError();

  bestHit_wrapper(stream, layer, d_XHitSize, d_XHitArr,
                  d_Err_iP, d_msErr, d_msPar, d_par_iP, d_outChi2,
                  d_Chi2, d_HitsIdx,
                  maxSize, N);
  kalmanUpdate_wrapper(stream, d_Err_iP, d_msErr,
                       d_par_iP, d_msPar, d_par_iC, d_Err_iC, N);
  if (ilay + 1 < Config::nLayers) {
    propagationForBuilding_wrapper(stream, d_Err_iC, d_par_iC, d_inChg, 
        radius, d_Err_iP, d_par_iP, N); 
  }
}   
#endif

#if 1
template <typename T>
void FitterCU<T>::propagateTracksToR(float radius, int N) {
  propagationForBuilding_wrapper(stream, d_Err_iC, d_par_iC, d_inChg, 
                                 radius, d_Err_iP, d_par_iP, N); 
}
#endif

template <typename T>
void FitterCU<T>::propagateTracksToR_standalone(float radius, int N,
    MPlexLS& Err_iC, MPlexLV& par_iC, MPlexQI& inChg, MPlexLS& Err_iP, MPlexLV& Par_iP) {
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

    d_msPar.copyAsyncFromHost(stream, msPar[hi]);
    d_msErr.copyAsyncFromHost(stream, msErr[hi]);

    propagationMerged();
    kalmanUpdateMerged();
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

