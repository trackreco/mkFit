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

#include <iostream>
#include <stdlib.h>
#include <sys/time.h>

// Set this to 8 for AVX, 16 for MIC
#ifdef __MIC__
const idx_t N = 16;
#else
const idx_t N = 8;
#endif

typedef ROOT::Math::SMatrix<float,6> SMatrix66;

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
  SMatrix66 res[N];
  SMatrix66 mul[N];
  SMatrix66 muz[N];
  for (int i=0; i<N; ++i)
    for (int j=0; j<6; ++j)
      for (int k=0; k<6; ++k)
      {
        res[i](j,k) = (float)rand() / RAND_MAX;
        mul[i](j,k) = (float)rand() / RAND_MAX;
        muz[i](j,k) = (float)rand() / RAND_MAX;
      }

  // ----------------------------------------------------------------

  double t0;
  t0 = dtime();

  for (int i = 0; i < 1000000; ++i)
  {
    //#pragma omp simd collapse(8)
    #pragma ivdep
    for (int m = 0; m < N; ++m)
    {
      res[m] = mul[m] * muz[m];
    }
  }

  double tsm = dtime() - t0;
  std::cout << "SMatrix run time = " << tsm << " s\n";

  // ----------------------------------------------------------------

  Matriplex<float, 6, 6, N> mpl, mpz, mpres;

  for (int i = 0; i < N; ++i)
  {
    mpl.Assign(i, mul[i].Array());
    mpz.Assign(i, muz[i].Array());
  }

  t0 = dtime();

  for (int i = 0; i < 1000000; ++i)
  {
    Multiply(mpl, mpz, mpres);
  }

  double tmp = dtime() - t0;
  std::cout << "Matriplex run time = " << tmp << " s\n";

  std::cout << "SMatrix / Matriplex = " << tsm/tmp << "\n";

  return 0;
}
