#include <cstdlib>

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
void FitterCU<T>::propagationMerged() {
  propagation_wrapper(stream, d_msPar, d_par_iC, d_inChg,
                      //d_par_iP, d_Err_iC, d_Err_iP, N); // TODO: Check outErr/errorProp
                      d_par_iP, d_errorProp, d_Err_iP, N);
}

template <typename T>
void FitterCU<T>::computeChi2gpu(
    const MPlexLS &psErr, MPlexHS &msErr,
    MPlexHV& msPar, const MPlexLV& propPar, GPlexQF& d_outChi2, int NN) {

  // TODO: add CMSGeom
  if (Config::useCMSGeom) {
    //propagateHelixToRMPlex(psErr,  psPar, inChg,  msPar, propErr, propPar);
    throw std::runtime_error("useCMSGeom not implemented yet for GPU");
  } else {}

  computeChi2_wrapper(stream, d_Err_iP, d_msErr, //d_resErr, 
      d_msPar, d_par_iP, d_outChi2, N);
}

// FIXME: Temporary. Separate allocations / transfers
template <typename T>
void FitterCU<T>::allocate_extra_addBestHit() {
  d_outChi2.allocate(Nalloc, QF);
  d_XHitPos.allocate(Nalloc, QI);
  d_XHitSize.allocate(Nalloc, QI);
  // FIXME: Make those GPlex-es. and use .allocate()
  cudaMalloc((void**)&d_HitsIdx, Nalloc*sizeof(int)); cudaCheckError();
  cudaMalloc((void**)&d_Chi2, Nalloc*sizeof(float)); cudaCheckError();
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::free_extra_addBestHit() {
  cudaFree(d_HitsIdx); cudaCheckError();
  cudaFree(d_Chi2); cudaCheckError();

  d_XHitPos.free(); cudaCheckError();
  d_XHitSize.free(); cudaCheckError();
  d_outChi2.free(); cudaCheckError();
}

// FIXME: Temporary. Separate allocations / transfers
template <typename T>
void FitterCU<T>::prepare_addBestHit(
    const MPlexLS &psErr, const MPlexLV& propPar,
    const MPlexQI &inChg, 
    size_t NN) {
  setNumberTracks(NN);  // temporary: should be end - beg

  createStream();
  cudaCheckError()

  // psErr -> d_Err_iP
  cudaMemcpy2DAsync(d_Err_iP.ptr, d_Err_iP.pitch, psErr.fArray, N*sizeof(T),
               N*sizeof(T), LS, cudaMemcpyHostToDevice, stream);
  // sendOutParToDevice(propPar);  // d_par_iP
  d_par_iP.copyAsyncFromHost(stream, propPar);
  //sendInChgToDevice(inChg);
  d_inChg.copyAsyncFromHost(stream, inChg);
}

// TODO: Temporary. Separate allocations / transfers
template <typename T>
void FitterCU<T>::finalize_addBestHit(
    MPlexHS &msErr, MPlexHV& msPar,
    MPlexLS &outErr, MPlexLV &outPar,
    MPlexQI &HitsIdx, MPlexQF &Chi2) {
  //getOutParFromDevice(outPar);  // <- d_par_iC
  d_par_iC.copyAsyncToHost(stream, outPar);
  //getOutErrFromDevice(outErr);  // <- d_Err_iC
  d_Err_iC.copyAsyncToHost(stream, outErr);

  //
  // Get msPar, msErr, chi2 and HitIdx out from the GPU to the CPU
  cudaMemcpy2DAsync(msPar.fArray, N*sizeof(T), d_msPar.ptr, d_msPar.pitch, 
               N*sizeof(T), HV, cudaMemcpyDeviceToHost, stream);
  cudaMemcpy2DAsync(msErr.fArray, N*sizeof(T), d_msErr.ptr, d_msErr.pitch, 
               N*sizeof(T), HS, cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(HitsIdx.fArray, d_HitsIdx, N*sizeof(int), cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(Chi2.fArray, d_Chi2, N*sizeof(float), cudaMemcpyDeviceToHost, stream);


  destroyStream();
}

template <typename T>
void FitterCU<T>::addBestHit(BunchOfHitsCU &bunch) {

  selectHitRanges_wrapper(stream, bunch, d_XHitPos, d_XHitSize, 
      d_Err_iP, d_par_iP, N);

  // TODO: get this thing inside bestHit_kernel
  int maxSize = getMaxNumHits_wrapper(d_XHitSize, N);

  bestHit_wrapper(stream, bunch, d_XHitPos,
                  d_Err_iP, d_msErr, d_msPar, d_par_iP, d_outChi2,
                  d_Chi2, d_HitsIdx,
                  maxSize, N);

  // updateParametersMPlex
  kalmanUpdate_wrapper(stream, d_Err_iP, d_msErr,
                       d_par_iP, d_msPar, d_par_iC, d_Err_iC, N);
  //updateParametersMPlex(Err[iP], Par[iP], Chg, msErr[Nhits], msPar[Nhits],
	//    Err[iC], Par[iC]);
}   

#if 0 
template <typename T>
void FitterCU<T>::propagateTracksToR(float radius, int N) {
  //propagateHelixToRMPlex(Err[iC], Par[iC], Chg, R,
                         //Err[iP], Par[iP], N_proc);
  propagationForBuilding_wrapper(stream, radius,
    d_par_iC, d_inChg, d_par_iP, d_errorProp, d_Err_iP, N);
  //propagation_wrapper(stream, d_msPar, d_par_iC, d_inChg,
  //                    d_par_iP, d_errorProp, d_Err_iP, N);
}
#endif

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

  //sendInChgToDevice(Chg);
  d_inChg.copyAsyncFromHost(stream, Chg);
  //sendInParToDevice(par_iC);
  d_par_iC.copyAsyncFromHost(stream, par_iC);
  //sendInErrToDevice(err_iC);
  d_Err_iC.copyAsyncFromHost(stream, err_iC);

  cudaEventRecord(start, 0);
 
  double total_reorg = 0.;
  for (int hi = 0; hi < Nhits; ++hi)
  {
    // Switch outPut and inPut parameters and errors
    // similar to iC <-> iP in the CPU code.
    //setOutParFromInPar();
    d_par_iP.copyAsyncFromDevice(stream, d_par_iC); 
    //setOutErrFromInErr(); // d_Err_iP
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

    //sendMsParToDevice(msPar[hi]);
    d_msPar.copyAsyncFromHost(stream, msPar[hi]);
    //sendMsErrToDevice(msErr[hi]);
    d_msErr.copyAsyncFromHost(stream, msErr[hi]);

    propagationMerged();
    kalmanUpdateMerged();
  }
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  cudaEventElapsedTime(&etime, start, stop);
  std::cerr << "CUDA etime: " << etime << " ms.\n";
  std::cerr << "Total reorg: " << total_reorg << " ms.\n";

  //getOutParFromDevice(par_iC);
  d_par_iC.copyAsyncToHost(stream, par_iC);
  //getOutErrFromDevice(err_iC);
  d_Err_iC.copyAsyncToHost(stream, err_iC);
  
  cudaStreamSynchronize(stream);
  // freeDevice(); -> moved to mkFit/mkFit.cc
  destroyStream();

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}

