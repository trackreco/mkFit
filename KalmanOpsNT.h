void MultForKalmanGain(const MatriplexSym<float, 6>& A,
                       const MatriplexSym<float, 6>& B,
                             Matriplex<float, 6, 6>& C)
{
   // calculate Kalman gain -- multiplication where B / resErr is only populated
   // in upper-left 3x3
   // kalmanGain = propErr * resErrInverse
   //     C      =    A    *      B

   const idx_t N = A.N;

  __assume_aligned(C.fArray, 64);
  __assume(N%16==0);

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      C.fArray[0 * N + n] = A.fArray[0 * N + n] * B.fArray[0 * N + n] + A.fArray[1 * N + n] * B.fArray[1 * N + n] + A.fArray[3 * N + n] * B.fArray[3 * N + n];
      C.fArray[1 * N + n] = A.fArray[0 * N + n] * B.fArray[1 * N + n] + A.fArray[1 * N + n] * B.fArray[2 * N + n] + A.fArray[3 * N + n] * B.fArray[4 * N + n];
      C.fArray[2 * N + n] = A.fArray[0 * N + n] * B.fArray[3 * N + n] + A.fArray[1 * N + n] * B.fArray[4 * N + n] + A.fArray[3 * N + n] * B.fArray[5 * N + n];
      C.fArray[6 * N + n] = A.fArray[1 * N + n] * B.fArray[0 * N + n] + A.fArray[2 * N + n] * B.fArray[1 * N + n] + A.fArray[4 * N + n] * B.fArray[3 * N + n];
      C.fArray[7 * N + n] = A.fArray[1 * N + n] * B.fArray[1 * N + n] + A.fArray[2 * N + n] * B.fArray[2 * N + n] + A.fArray[4 * N + n] * B.fArray[4 * N + n];
      C.fArray[8 * N + n] = A.fArray[1 * N + n] * B.fArray[3 * N + n] + A.fArray[2 * N + n] * B.fArray[4 * N + n] + A.fArray[4 * N + n] * B.fArray[5 * N + n];
      C.fArray[12 * N + n] = A.fArray[3 * N + n] * B.fArray[0 * N + n] + A.fArray[4 * N + n] * B.fArray[1 * N + n] + A.fArray[5 * N + n] * B.fArray[3 * N + n];
      C.fArray[13 * N + n] = A.fArray[3 * N + n] * B.fArray[1 * N + n] + A.fArray[4 * N + n] * B.fArray[2 * N + n] + A.fArray[5 * N + n] * B.fArray[4 * N + n];
      C.fArray[14 * N + n] = A.fArray[3 * N + n] * B.fArray[3 * N + n] + A.fArray[4 * N + n] * B.fArray[4 * N + n] + A.fArray[5 * N + n] * B.fArray[5 * N + n];
      C.fArray[18 * N + n] = A.fArray[6 * N + n] * B.fArray[0 * N + n] + A.fArray[7 * N + n] * B.fArray[1 * N + n] + A.fArray[8 * N + n] * B.fArray[3 * N + n];
      C.fArray[19 * N + n] = A.fArray[6 * N + n] * B.fArray[1 * N + n] + A.fArray[7 * N + n] * B.fArray[2 * N + n] + A.fArray[8 * N + n] * B.fArray[4 * N + n];
      C.fArray[20 * N + n] = A.fArray[6 * N + n] * B.fArray[3 * N + n] + A.fArray[7 * N + n] * B.fArray[4 * N + n] + A.fArray[8 * N + n] * B.fArray[5 * N + n];
      C.fArray[24 * N + n] = A.fArray[10 * N + n] * B.fArray[0 * N + n] + A.fArray[11 * N + n] * B.fArray[1 * N + n] + A.fArray[12 * N + n] * B.fArray[3 * N + n];
      C.fArray[25 * N + n] = A.fArray[10 * N + n] * B.fArray[1 * N + n] + A.fArray[11 * N + n] * B.fArray[2 * N + n] + A.fArray[12 * N + n] * B.fArray[4 * N + n];
      C.fArray[26 * N + n] = A.fArray[10 * N + n] * B.fArray[3 * N + n] + A.fArray[11 * N + n] * B.fArray[4 * N + n] + A.fArray[12 * N + n] * B.fArray[5 * N + n];
      C.fArray[30 * N + n] = A.fArray[15 * N + n] * B.fArray[0 * N + n] + A.fArray[16 * N + n] * B.fArray[1 * N + n] + A.fArray[17 * N + n] * B.fArray[3 * N + n];
      C.fArray[31 * N + n] = A.fArray[15 * N + n] * B.fArray[1 * N + n] + A.fArray[16 * N + n] * B.fArray[2 * N + n] + A.fArray[17 * N + n] * B.fArray[4 * N + n];
      C.fArray[32 * N + n] = A.fArray[15 * N + n] * B.fArray[3 * N + n] + A.fArray[16 * N + n] * B.fArray[4 * N + n] + A.fArray[17 * N + n] * B.fArray[5 * N + n];
   }
}

//------------------------------------------------------------------------------

void MultResidualsAdd(const Matriplex<float, 6, 6>& A,
                      const Matriplex<float, 6, 1>& B,
                      const Matriplex<float, 6, 1>& C,
                            Matriplex<float, 6, 1>& D)
{
   // outPar = psPar + kalmanGain*(msPar-psPar)
   //   D    =   B         A         C  -  B
   // where right half of kalman gain is 0 

   const idx_t N = A.N;

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      // manually subrtact into local vars -- 3 of them
      float d0 = C.fArray[0 * N + n] - B.fArray[0 * N + n];
      float d1 = C.fArray[1 * N + n] - B.fArray[1 * N + n];
      float d2 = C.fArray[2 * N + n] - B.fArray[2 * N + n];

      // generate loop (can also write it manually this time, it's not much)
      D.fArray[0 * N + n] = B.fArray[0 * N + n] + A.fArray[ 0 * N + n] * d0 + A.fArray[ 1 * N + n] * d1 + A.fArray[ 2 * N + n] * d2;
      D.fArray[1 * N + n] = B.fArray[1 * N + n] + A.fArray[ 6 * N + n] * d0 + A.fArray[ 7 * N + n] * d1 + A.fArray[ 8 * N + n] * d2;
      D.fArray[2 * N + n] = B.fArray[2 * N + n] + A.fArray[12 * N + n] * d0 + A.fArray[13 * N + n] * d1 + A.fArray[14 * N + n] * d2;
      D.fArray[3 * N + n] = B.fArray[3 * N + n] + A.fArray[18 * N + n] * d0 + A.fArray[19 * N + n] * d1 + A.fArray[20 * N + n] * d2;
      D.fArray[4 * N + n] = B.fArray[4 * N + n] + A.fArray[24 * N + n] * d0 + A.fArray[25 * N + n] * d1 + A.fArray[26 * N + n] * d2;
      D.fArray[5 * N + n] = B.fArray[5 * N + n] + A.fArray[30 * N + n] * d0 + A.fArray[31 * N + n] * d1 + A.fArray[32 * N + n] * d2;
   }
}

//------------------------------------------------------------------------------

void FinalKalmanErr(const MatriplexSym<float, 6>& A,
                    const Matriplex<float, 6, 6>& B,
                          MatriplexSym<float, 6>& C)
{
   // propErr - kalmanGain*propErr
   //    A    -      B    *   A
   // where right half of B is 0

   const idx_t N = A.N;

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      C.fArray[0 * N + n] = A.fArray[0 * N + n] - B.fArray[0 * N + n] * A.fArray[0 * N + n] - B.fArray[1 * N + n] * A.fArray[1 * N + n] - B.fArray[2 * N + n] * A.fArray[3 * N + n];
      C.fArray[1 * N + n] = A.fArray[1 * N + n] - B.fArray[6 * N + n] * A.fArray[0 * N + n] - B.fArray[7 * N + n] * A.fArray[1 * N + n] - B.fArray[8 * N + n] * A.fArray[3 * N + n];
      C.fArray[2 * N + n] = A.fArray[2 * N + n] - B.fArray[6 * N + n] * A.fArray[1 * N + n] - B.fArray[7 * N + n] * A.fArray[2 * N + n] - B.fArray[8 * N + n] * A.fArray[4 * N + n];
      C.fArray[3 * N + n] = A.fArray[3 * N + n] - B.fArray[12 * N + n] * A.fArray[0 * N + n] - B.fArray[13 * N + n] * A.fArray[1 * N + n] - B.fArray[14 * N + n] * A.fArray[3 * N + n];
      C.fArray[4 * N + n] = A.fArray[4 * N + n] - B.fArray[12 * N + n] * A.fArray[1 * N + n] - B.fArray[13 * N + n] * A.fArray[2 * N + n] - B.fArray[14 * N + n] * A.fArray[4 * N + n];
      C.fArray[5 * N + n] = A.fArray[5 * N + n] - B.fArray[12 * N + n] * A.fArray[3 * N + n] - B.fArray[13 * N + n] * A.fArray[4 * N + n] - B.fArray[14 * N + n] * A.fArray[5 * N + n];
      C.fArray[6 * N + n] = A.fArray[6 * N + n] - B.fArray[18 * N + n] * A.fArray[0 * N + n] - B.fArray[19 * N + n] * A.fArray[1 * N + n] - B.fArray[20 * N + n] * A.fArray[3 * N + n];
      C.fArray[7 * N + n] = A.fArray[7 * N + n] - B.fArray[18 * N + n] * A.fArray[1 * N + n] - B.fArray[19 * N + n] * A.fArray[2 * N + n] - B.fArray[20 * N + n] * A.fArray[4 * N + n];
      C.fArray[8 * N + n] = A.fArray[8 * N + n] - B.fArray[18 * N + n] * A.fArray[3 * N + n] - B.fArray[19 * N + n] * A.fArray[4 * N + n] - B.fArray[20 * N + n] * A.fArray[5 * N + n];
      C.fArray[9 * N + n] = A.fArray[9 * N + n] - B.fArray[18 * N + n] * A.fArray[6 * N + n] - B.fArray[19 * N + n] * A.fArray[7 * N + n] - B.fArray[20 * N + n] * A.fArray[8 * N + n];
      C.fArray[10 * N + n] = A.fArray[10 * N + n] - B.fArray[24 * N + n] * A.fArray[0 * N + n] - B.fArray[25 * N + n] * A.fArray[1 * N + n] - B.fArray[26 * N + n] * A.fArray[3 * N + n];
      C.fArray[11 * N + n] = A.fArray[11 * N + n] - B.fArray[24 * N + n] * A.fArray[1 * N + n] - B.fArray[25 * N + n] * A.fArray[2 * N + n] - B.fArray[26 * N + n] * A.fArray[4 * N + n];
      C.fArray[12 * N + n] = A.fArray[12 * N + n] - B.fArray[24 * N + n] * A.fArray[3 * N + n] - B.fArray[25 * N + n] * A.fArray[4 * N + n] - B.fArray[26 * N + n] * A.fArray[5 * N + n];
      C.fArray[13 * N + n] = A.fArray[13 * N + n] - B.fArray[24 * N + n] * A.fArray[6 * N + n] - B.fArray[25 * N + n] * A.fArray[7 * N + n] - B.fArray[26 * N + n] * A.fArray[8 * N + n];
      C.fArray[14 * N + n] = A.fArray[14 * N + n] - B.fArray[24 * N + n] * A.fArray[10 * N + n] - B.fArray[25 * N + n] * A.fArray[11 * N + n] - B.fArray[26 * N + n] * A.fArray[12 * N + n];
      C.fArray[15 * N + n] = A.fArray[15 * N + n] - B.fArray[30 * N + n] * A.fArray[0 * N + n] - B.fArray[31 * N + n] * A.fArray[1 * N + n] - B.fArray[32 * N + n] * A.fArray[3 * N + n];
      C.fArray[16 * N + n] = A.fArray[16 * N + n] - B.fArray[30 * N + n] * A.fArray[1 * N + n] - B.fArray[31 * N + n] * A.fArray[2 * N + n] - B.fArray[32 * N + n] * A.fArray[4 * N + n];
      C.fArray[17 * N + n] = A.fArray[17 * N + n] - B.fArray[30 * N + n] * A.fArray[3 * N + n] - B.fArray[31 * N + n] * A.fArray[4 * N + n] - B.fArray[32 * N + n] * A.fArray[5 * N + n];
      C.fArray[18 * N + n] = A.fArray[18 * N + n] - B.fArray[30 * N + n] * A.fArray[6 * N + n] - B.fArray[31 * N + n] * A.fArray[7 * N + n] - B.fArray[32 * N + n] * A.fArray[8 * N + n];
      C.fArray[19 * N + n] = A.fArray[19 * N + n] - B.fArray[30 * N + n] * A.fArray[10 * N + n] - B.fArray[31 * N + n] * A.fArray[11 * N + n] - B.fArray[32 * N + n] * A.fArray[12 * N + n];
      C.fArray[20 * N + n] = A.fArray[20 * N + n] - B.fArray[30 * N + n] * A.fArray[15 * N + n] - B.fArray[31 * N + n] * A.fArray[16 * N + n] - B.fArray[32 * N + n] * A.fArray[17 * N + n];
   }
}

//------------------------------------------------------------------------------

