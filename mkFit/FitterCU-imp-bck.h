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
  d_par_iC.allocate(N, LV);
  d_inChg.allocate(N, QI);
  d_par_iP.allocate(N, LV);
  d_msRad.allocate(N, QF);
  d_errorProp.allocate(N, LL);
  d_Err_iP.allocate(N, LS);
  d_msPar.allocate(N, HV);
  d_kalmanGain.allocate(N, LH);
  d_outErr.allocate(N, LS);
  d_msErr.allocate(N, HS);
  d_resErr.allocate(N, HS);

  cudaCheckError()
}

template <typename T>
void FitterCU<T>::freeDevice() {
  d_par_iC.free();
  d_inChg.free();
  d_par_iP.free();
  d_msRad.free();
  d_errorProp.free();
  d_Err_iP.free();
  d_msPar.free();
  d_kalmanGain.free();
  d_outErr.free();
  d_msErr.free();
  d_resErr.free();

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
void FitterCU<T>::swap_iP_iC() {
  T *tmp;
  tmp = d_par_iC.ptr;
  d_par_iC.ptr = d_par_iP.ptr;
  d_par_iP.ptr = tmp;
}

template <typename T>
void FitterCU<T>::computeMsRad() {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  computeMsRad_wrapper(grid, block, stream, d_msPar, d_msRad, N);
 
  cudaCheckErrorSync();
}

template <typename T>
void FitterCU<T>::helixAtRFromIterative() {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  helixAtRFromIterative_wrapper(grid, block, stream, d_par_iC, d_inChg, d_par_iP,
      d_msRad, d_errorProp, N);
  cudaCheckErrorSync();
}

template <typename T>
void FitterCU<T>::similarity() {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  similarity_wrapper(grid, block, stream, d_errorProp, d_Err_iP, N);
  cudaCheckErrorSync();
}

///////////////////////////////////////////////////////////////////////////////
///  Updater-specific methods                                               ///
///////////////////////////////////////////////////////////////////////////////

template <typename T>
void FitterCU<T>::multKalmanGainCU() {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  upParam_MultKalmanGain_wrapper(grid, block, stream,
      d_Err_iP, d_resErr, d_kalmanGain, N);
  cudaCheckErrorSync();
}


template <typename T>
void FitterCU<T>::multResidualsAdd() {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  
  multResidualsAdd_wrapper(grid, block, stream,
      d_kalmanGain, d_par_iP, d_msPar, d_par_iC, N);
  cudaCheckErrorSync();
}


template <typename T>
void FitterCU<T>::kalmanGain_x_propErr() {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  
  kalmanGain_x_propErr_wrapper(grid, block, stream,
      d_kalmanGain, d_Err_iP, d_outErr, N);
  cudaCheckErrorSync();
}


template <typename T>
void FitterCU<T>::InvertCramerSym() {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  invertCramerSym_wrapper(grid, block, stream, d_resErr, N);
  cudaCheckErrorSync();
}

template <typename T>
void FitterCU<T>::addIntoUpperLeft3x3() {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);

  addIntoUpperLeft3x3_wrapper(grid, block, stream, d_Err_iP, d_msErr, d_resErr, N);
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

  createStream();

  allocateDevice();
  sendInErrToDevice(err_iC);
  sendInParToDevice(par_iC);
  sendInChgToDevice(Chg);


  cudaEventRecord(start, 0);
#if 0
  // try to reorganize on the gpu
  float *inErr = new float[LS*N];
  float *inPar = new float[LV*N];
  int *inChg = new float[QI*N];
  for (int i = beg; i < end; ++i) {
    memcpy(inErr + i*LS, tracks[i].errors().Array(), LS*sizeof(float));
    memcpy(inPar + i*LV, tracks[i].parameters().Array(), LV*sizeof(float));
    inChg[i*QI] = tracks[i].charge();
  }
  sendInChgToDevice(inChg);
#endif

#if 0
  GPlex<float> t_inErr;
  GPlex<float> t_inPar;
  t_inErr.allocate(LS, N);
  t_inPar.allocate(LV, N);
  cudaCheckError();
  int itrack = 0;
  cudaMemcpy2DAsync(t_inErr.ptr, t_inErr.pitch,
                    inErr, LS*sizeof(float),
                    LS*sizeof(float), N,
                    cudaMemcpyHostToDevice, stream);
  cudaMemcpy2DAsync(t_inPar.ptr, t_inPar.pitch,
                    inPar, LV*sizeof(float),
                    LV*sizeof(float), N,
                    cudaMemcpyHostToDevice, stream);

  toMatriplex_wrapper(stream, d_outErr, t_inErr, N, LS);
  toMatriplex_wrapper(stream, d_par_iC, t_inPar, N, LV);
  //cudaDeviceSynchronize();
  //exit(-1);
#endif
 
  for (int hi = 0; hi < Nhits; ++hi)
  {
    setOutParFromInPar();
    setOutErrFromInErr(); // d_Err_iP

#if 0
    int itrack = 0;
    for (int i = beg; i < end; ++i, ++itrack)
    {
      Track &trk = tracks[i];

      const int hidx = trk.getHitIdx(hi);
      const Hit &hit = layerHits[hi][hidx];

      msErr[hi].CopyIn(itrack, hit.errArray());
      msPar[hi].CopyIn(itrack, hit.posArray());
    }
    //cudaDeviceSynchronize(); 
#endif

    sendMsParToDevice(msPar[hi]);
    sendMsErrToDevice(msErr[hi]);

    computeMsRad();
    if (Config::doIterative) {
      helixAtRFromIterative();
    } else {
      ; // Nothing for now
    }

    similarity();

    addIntoUpperLeft3x3();
    InvertCramerSym();
    multKalmanGainCU();
    multResidualsAdd();
    kalmanGain_x_propErr();
  }
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  cudaEventElapsedTime(&etime, start, stop);
  std::cerr << "CUDA etime: " << etime << " ms.\n";

  getOutParFromDevice(par_iC);
  getOutErrFromDevice(err_iC);
  
  cudaDeviceSynchronize();  // TODO: change with cudaStreamSynchronize();
  freeDevice();

  //t_inErr.free();
  //t_inPar.free();
  //delete[] inErr;
  //delete[] inPar;
  //delete[] inChg;

  cudaCheckError();
  destroyStream();

  cudaEventDestroy(start);
  cudaEventDestroy(stop);

}
