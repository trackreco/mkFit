// HOST / AVX
// icc -O3 -mavx -I. -o mult66 mult66.cc  -vec-report=3
// ./mult66

// MIC
// icc -O3 -mmic -I. -o mult66-mic mult66.cc -vec-report=3
// scp mult66-mic root@mic0:
// ssh root@mic0 ./mult66-mic

#include "Matriplex.h"

#define ROOT_Math_MnConfig
#include "Math/SMatrix.h"

#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <sys/time.h>

// Set this to 8 for AVX, 16 for MIC
#ifdef __MIC__
const idx_t N = 16;
#else
const idx_t N = 8;
#endif

#ifdef MDIM
const idx_t M = MDIM;
#else
const idx_t M = 6;
#endif

const int NN_MULT = 1000000;
const int NN_INV  = 1;//000000;

typedef ROOT::Math::SMatrix<float, M> SMatrixMM;

double dtime()
{
    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday(&mytime,(struct timezone*)0);
    tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
    return( tseconds );
}

int main()
{
  srand(23545);
  SMatrixMM res[N];
  SMatrixMM mul[N];
  SMatrixMM muz[N];
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < M; ++j)
      for (int k = 0; k < M; ++k)
      {
        mul[i](j,k) = (float)rand() / RAND_MAX;
        muz[i](j,k) = (float)rand() / RAND_MAX;
      }

  // ================================================================

  double t0;
  t0 = dtime();

  for (int i = 0; i < NN_MULT; ++i)
  {
    //#pragma omp simd collapse(8)
#pragma ivdep
    for (int m = 0; m < N; ++m)
    {
      res[m] = mul[m] * muz[m];
    }
  }

  double tsm = dtime() - t0;
  std::cout << "SMatrix multiply time = " << tsm << " s\n";

  // ----------------------------------------------------------------

  Matriplex<float, M, M, N> mpl, mpz, mpres;

  for (int i = 0; i < N; ++i)
  {
    mpl.Assign(i, mul[i].Array());
    mpz.Assign(i, muz[i].Array());
  }

  t0 = dtime();

  for (int i = 0; i < NN_MULT; ++i)
  {
    Multiply(mpl, mpz, mpres);
  }

  double tmp = dtime() - t0;
  std::cout << "Matriplex multiply time = " << tmp << " s\n";

  std::cout << "SMatrix / Matriplex = " << tsm/tmp << "\n\n";

  for (int i = 0; i < N; ++i)
    for (int j = 0; j < M; ++j)
      for (int k = 0; k < M; ++k)
      {
        if (res[i](j,k) != mpres.At(j, k, i))
          std::cout << i <<" "<< j <<" "<< k <<" "<< res[i](j,k) <<" "<< mpres.At(j, k, i) << "\n";
      }

  // ================================================================

  t0 = dtime();

  for (int i = 0; i < NN_INV; ++i)
  {
    //#pragma omp simd collapse(8)
#pragma ivdep
    for (int m = 0; m < N; ++m)
    {
      mul[m].InvertFast();
      muz[m].InvertFast();
    }
  }
  tsm = dtime() - t0;
  std::cout << "SMatrix invert time = " << tsm << " s\n";

  // ----------------------------------------------------------------

  t0 = dtime();

  for (int i = 0; i < NN_INV; ++i)
  {
    InvertCramer(mpl);
    InvertCramer(mpz);
  }

  tmp = dtime() - t0;
  std::cout << "Matriplex invert time = " << tmp << " s\n";
  std::cout << "SMatrix / Matriplex = " << tsm/tmp << "\n\n";

  // ================================================================
  for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
      printf("%8f ", mul[0](i,j)); printf("\n");
  } printf("\n");

  for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
      printf("%8f ", mpl.At(i,j,0)); printf("\n");
  } printf("\n");

  return 0;
}



#ifdef ZMAJLA_ZAJLA

  std::cout << "LLLLL\n";

  for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
      printf("%8f ", mul[0](i,j)); printf("\n");
  } printf("\n");
  for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
      printf("%8f ", mpl.At(i,j,0)); printf("\n");
  } printf("\n");
  for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
      printf("%8f ", muz[0](i,j)); printf("\n");
  } printf("\n");
  for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
      printf("%8f ", mpz.At(i,j,0)); printf("\n");
  } printf("\n");

  for (int i = 0; i < N; ++i)
    for (int j = 0; j < M; ++j)
      for (int k = 0; k < M; ++k)
      {
        if (mul[i](j,k) != mpl.At(j, k, i))
          std::cout << i <<" "<< j <<" "<< k <<" "<< mul[i](j,k) <<" "<< mpl.At(j, k, i) << "\n";
      }

  std::cout << "ZZZZZ\n";

  for (int i = 0; i < N; ++i)
    for (int j = 0; j < M; ++j)
      for (int k = 0; k < M; ++k)
      {
        if (muz[i](j,k) != mpz.At(j, k, i))
          std::cout << i <<" "<< j <<" "<< k <<" "<< muz[i](j,k) <<" "<< mpz.At(j, k, i) << "\n";
      }

  std::cout << "RESRES\n";

#endif
