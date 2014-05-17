#define ROOT_Math_MnConfig
#include "Math/SMatrix.h"

#include <random>

#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <sys/time.h>

// Set this to 8 for AVX, 16 for MIC
const idx_t Sfac = 1;
#ifdef __MIC__
const idx_t S = 16 * Sfac;
#else
const idx_t S = 8  * Sfac;
#endif

#ifdef MDIM
const idx_t M = MDIM;
#else
const idx_t M = 6;
#endif

const int N = 256 * 10;

const int Nm = N / S;
const int Ns = N;

const int NN_MULT = 10000;
const int NN_INV  = 0;

const int ALIGN   = 4096;

const int COMPARE_INVERSES = 0;

#define SYMM 1

typedef ROOT::Math::SMatrix<float, M>                                      SMatrixMM;
typedef ROOT::Math::SMatrix<float, M, M, ROOT::Math::MatRepSym<float, M> > SMatrixSS;



double dtime()
{
    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday(&mytime,(struct timezone*)0);
    tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
    return( tseconds );
}

template <typename X>
X* new_sth(int n)
{
  return (X*) _mm_malloc(sizeof(X) * n, ALIGN);
}

void init_mulmuz(SMatrixMM *mul,  SMatrixMM *muz,
                 SMatrixSS *muls, SMatrixSS *muzs)
{
  std::default_random_engine gen(0xbeef0133);
  std::normal_distribution<float> dis(5.0, 0.05);

  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < M; ++j)
    {
      for (int k = 0; k < M; ++k)
      {
        mul[i](j,k) = dis(gen);
        muz[i](j,k) = dis(gen);
      }
    }
    mul[i] = Transpose(mul[i]) * mul[i];
    muz[i] = Transpose(muz[i]) * muz[i];

    if (muls && muzs)
    {
      for (int j = 0; j < M; ++j)
      {
        for (int k = j; k < M; ++k)
        {
          muls[i](j,k) = mul[i](j,k);
          muzs[i](j,k) = muz[i](j,k);
        }
      }
    }
  }
}

void init_mulmuz(SMatrixSS *muls, SMatrixSS *muzs)
{
  std::default_random_engine gen(0xbeef0133);
  std::normal_distribution<float> dis(5.0, 0.05);

  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < M; ++j)
    {
      for (int k = j; k < M; ++k)
      {
        muls[i](j,k) = dis(gen);
        muzs[i](j,k) = dis(gen);
      }
    }
  }
}
