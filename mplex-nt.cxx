// HOST / AVX
// icc -O3 -mavx -I. -o mult66 mult66.cc  -vec-report=3
// ./mult66

// MIC
// icc -O3 -mmic -I. -o mult66-mic mult66.cc -vec-report=3
// scp mult66-mic root@mic0:
// ssh root@mic0 ./mult66-mic

#include "MatriplexNT.h"

#include "mplex-common.h"


typedef Matriplex<float, M, M> MPlexMM;


int main(int arg, char *argv[])
{
  SMatrixMM *res  = new_sth<SMatrixMM>(Ns);
  SMatrixMM *mul  = new_sth<SMatrixMM>(Ns);
  SMatrixMM *muz  = new_sth<SMatrixMM>(Ns);

#ifdef SYMM
  SMatrixSS *muls = new_sth<SMatrixSS>(Ns);
  SMatrixSS *muzs = new_sth<SMatrixSS>(Ns);
#else
  SMatrixSS *muls = 0;
  SMatrixSS *muzs = 0;
#endif

  MPlexMM mpl(N);
  MPlexMM mpz(N);
  MPlexMM mpres(N);
  

  init_mulmuz(mul, muz, muls, muzs);

  // Multiplex arrays

  for (int i = 0; i < N; ++i)
  {
    mpl.Assign(i, mul[i].Array());
    mpz.Assign(i, muz[i].Array());
  }

  // ================================================================

  double t0, tmp, tsm;

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

  tsm = dtime() - t0;
  std::cout << "SMatrix multiply time = " << tsm << " s\n";

  // ----------------------------------------------------------------

  t0 = dtime();

  for (int i = 0; i < NN_MULT; ++i)
  {
    Multiply(mpl, mpz, mpres);
  }


  tmp = dtime() - t0;
  std::cout << "Matriplex multiply time = " << tmp << " s\n";

  std::cout << "SMatrix / Matriplex = " << tsm/tmp << "\n\n";

  /*
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < M; ++j)
      for (int k = 0; k < M; ++k)
      {
    if (res[i](j,k) != mpres.At(j, k, i))
      std::cout << i <<" "<< j <<" "<< k <<" "<< res[i](j,k) <<" "<< mpres.At(j, k, i) << "\n";
  }
  */

  // ================================================================

  t0 = dtime();

  for (int i = 0; i < NN_INV; ++i)
  {
    //#pragma omp simd collapse(8)
#pragma ivdep
    for (int m = 0; m < N; ++m)
    {
    //mul[m].InvertFast();
    //muz[m].InvertFast();
#ifdef SYMM
    bool bl = muls[m].InvertChol();
    bool bz = muzs[m].InvertChol();
#else
    bool bl = mul[m].InvertChol();
    bool bz = muz[m].InvertChol();
#endif
    //if ( ! bl || ! bz)   printf("Grr %d %d %d\n", m, bl, bz);
    }
  }
  tsm = dtime() - t0;
  std::cout << "SMatrix invert time = " << tsm << " s\n";

  // ----------------------------------------------------------------

  t0 = dtime();

  for (int i = 0; i < NN_INV; ++i)
  {
    // InvertCramer(mpl);
    // InvertCramer(mpz);
    InvertChol(mpl);
    InvertChol(mpz);
  }

  tmp = dtime() - t0;
  std::cout << "Matriplex invert time = " << tmp << " s\n";
  std::cout << "SMatrix / Matriplex = " << tsm/tmp << "\n\n";

  // ================================================================
  for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
      printf("%8f ", muls[0](i,j)); printf("\n");
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
