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
  d_par_iC.allocate(Nalloc, LV);
  d_inChg.allocate(Nalloc, QI);
  d_par_iP.allocate(Nalloc, LV);
  d_errorProp.allocate(Nalloc, LL);
  d_Err_iP.allocate(Nalloc, LS);
  d_msPar.allocate(Nalloc, HV);
  d_outErr.allocate(Nalloc, LS);
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
  d_outErr.free();
  d_msErr.free();

  cudaCheckError()
}

template <typename T>
void FitterCU<T>::sendInParToDevice(const MPlexLV& inPar) {
  cudaMemcpy2DAsync(d_par_iC.ptr, d_par_iC.pitch, inPar.fArray, N*sizeof(T),
               N*sizeof(T), LV, cudaMemcpyHostToDevice, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::sendInErrToDevice(const MPlexLS& inErr) {
  cudaMemcpy2DAsync(d_outErr.ptr, d_outErr.pitch, inErr.fArray, N*sizeof(T),
               N*sizeof(T), LS, cudaMemcpyHostToDevice, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::sendInChgToDevice(const MPlexQI& inChg) {
  cudaMemcpy2DAsync(d_inChg.ptr, d_inChg.pitch, inChg.fArray, N*sizeof(T),
               N*sizeof(T), QI, cudaMemcpyHostToDevice, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::sendMsRadToDevice(const MPlexQF& msRad) {
  cudaMemcpy2DAsync(d_msRad.ptr, d_msRad.pitch, msRad.fArray, N*sizeof(T),
               N*sizeof(T), QF, cudaMemcpyHostToDevice, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::sendOutParToDevice(const MPlexLV& outPar) {
  cudaMemcpy2DAsync(d_par_iP.ptr, d_par_iP.pitch, outPar.fArray, N*sizeof(T),
               N*sizeof(T), LV, cudaMemcpyHostToDevice, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::sendOutErrToDevice(const MPlexLS& outErr) {
  cudaMemcpy2DAsync(d_Err_iP.ptr, d_Err_iP.pitch, outErr.fArray, N*sizeof(T),
               N*sizeof(T), LS, cudaMemcpyHostToDevice, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::sendMsParToDevice(const MPlexHV& msPar) {
  cudaMemcpy2DAsync(d_msPar.ptr, d_msPar.pitch, msPar.fArray, N*sizeof(T),
               N*sizeof(T), HV, cudaMemcpyHostToDevice, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::sendMsErrToDevice(const MPlexHS& msErr) {
  cudaMemcpy2DAsync(d_msErr.ptr, d_msErr.pitch, msErr.fArray, N*sizeof(T),
               N*sizeof(T), HS, cudaMemcpyHostToDevice, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::getOutParFromDevice(MPlexLV& outPar) {
  cudaMemcpy2DAsync(outPar.fArray, N*sizeof(T), d_par_iC.ptr, d_par_iC.pitch,
               N*sizeof(T), LV, cudaMemcpyDeviceToHost, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::getErrorPropFromDevice(MPlexLL& errorProp) {
  cudaMemcpy2DAsync(errorProp.fArray, N*sizeof(T),
               d_errorProp.ptr, d_errorProp.pitch,
               N*sizeof(T), LL, cudaMemcpyDeviceToHost, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::getOutErrFromDevice(MPlexLS& outErr) {
  cudaMemcpy2DAsync(outErr.fArray, N*sizeof(T), d_outErr.ptr, d_outErr.pitch,
               N*sizeof(T), LS, cudaMemcpyDeviceToHost, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::getMsRadFromDevice(MPlexQF& msRad) {
  cudaMemcpy2DAsync(msRad.fArray, N*sizeof(T), d_msRad.ptr, d_msRad.pitch,
               N*sizeof(T), QF, cudaMemcpyDeviceToHost, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::setOutParFromInPar() {
  cudaMemcpy2DAsync(d_par_iP.ptr, d_par_iP.pitch, d_par_iC.ptr, d_par_iC.pitch,
               N*sizeof(T), LV, cudaMemcpyDeviceToDevice, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::setOutErrFromInErr() {
  cudaMemcpy2DAsync(d_Err_iP.ptr, d_Err_iP.pitch, d_outErr.ptr, d_outErr.pitch,
               N*sizeof(T), LS, cudaMemcpyDeviceToDevice, stream);
  cudaCheckError()
}

template <typename T>
void FitterCU<T>::kalmanUpdateMerged() {
  kalmanUpdate_wrapper(stream, d_Err_iP, d_msErr,
                       d_par_iP, d_msPar, d_par_iC, d_outErr, N);
}

template <typename T>
void FitterCU<T>::propagationMerged() {
  propagation_wrapper(stream, d_msPar, d_par_iC, d_inChg,
                      d_par_iP, d_errorProp, d_Err_iP, N);
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

  sendInChgToDevice(Chg);
  sendInParToDevice(par_iC);
  sendInErrToDevice(err_iC);

  cudaEventRecord(start, 0);
 
  double total_reorg = 0.;
  for (int hi = 0; hi < Nhits; ++hi)
  {
    // Switch outPut and inPut parameters and errors
    // similar to iC <-> iP in the CPU code.
    setOutParFromInPar();
    setOutErrFromInErr(); // d_Err_iP
    
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

    sendMsParToDevice(msPar[hi]);
    sendMsErrToDevice(msErr[hi]);

    propagationMerged();
    kalmanUpdateMerged();
  }
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  cudaEventElapsedTime(&etime, start, stop);
  std::cerr << "CUDA etime: " << etime << " ms.\n";
  std::cerr << "Total reorg: " << total_reorg << " ms.\n";

  getOutParFromDevice(par_iC);
  getOutErrFromDevice(err_iC);
  
  cudaStreamSynchronize(stream);
  // freeDevice(); -> moved to mkFit/mkFit.cc
  destroyStream();

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}

