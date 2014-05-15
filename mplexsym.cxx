// HOST / AVX
// icc -O3 -mavx -I. -o mult66 mult66.cc  -vec-report=3
// ./mult66

// MIC
// icc -O3 -mmic -I. -o mult66-mic mult66.cc -vec-report=3
// scp mult66-mic root@mic0:
// ssh root@mic0 ./mult66-mic

#include "MatriplexSym.h"

#include "mplex-common.h"


typedef Matriplex<float, M, M, S> MPlexMM;
typedef MatriplexSym<float, M, S> MPlexSS;


int main()
{
  SMatrixSS *mul  = new_sth<SMatrixSS>(Ns);
  SMatrixSS *muz  = new_sth<SMatrixSS>(Ns);
  SMatrixMM *res  = new_sth<SMatrixMM>(Ns);

  MPlexSS *mpl   = new_sth<MPlexSS>(Nm);
  MPlexSS *mpz   = new_sth<MPlexSS>(Nm);
  MPlexMM *mpres = new_sth<MPlexMM>(Nm);


  init_mulmuz(mul, muz);

  // Multiplex arrays

  for (int i = 0; i < N; ++i)
  {
    mpl[i/S].Assign(i%S, mul[i].Array());
    mpz[i/S].Assign(i%S, muz[i].Array());
  }

  double t0, tmp, tsm;

  // ================================================================

  if (NN_MULT)
  {
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
      for (int j = 0; j < Nm; ++j)
      {
        Multiply(mpl[j], mpz[j], mpres[j]);
      }
    }


    tmp = dtime() - t0;
    std::cout << "Matriplex multiply time = " << tmp << " s\n";

    std::cout << "SMatrix / Matriplex = " << tsm/tmp << "";

    double x = 0, y = 0;
    for (int j = 0; j < Nm; ++j)
    {
      for (int k = 0; k < S; ++k)
      {
        x += res[j*S + k](1,2);
        y += mpres[j](1,2,k);
      }
    }
    std::cout << "\t\t\tx = " << x << ", y = " << y << "\n";

    std::cout << "\n";
    /*
      for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
      for (int k = 0; k < M; ++k)
      {
      if (res[i](j,k) != mpres.At(j, k, i))
      std::cout << i <<" "<< j <<" "<< k <<" "<< res[i](j,k) <<" "<< mpres.At(j, k, i) << "\n";
      }
    */
  }

  // ================================================================

  if (NN_INV)
  {
    t0 = dtime();

    for (int i = 0; i < NN_INV; ++i)
    {
      //#pragma omp simd collapse(8)
#pragma ivdep
      for (int m = 0; m < N; ++m)
      {
        //mul[m].InvertFast();
        //muz[m].InvertFast();
        bool bl = mul[m].InvertChol();
        bool bz = muz[m].InvertChol();

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
      for (int j = 0; j < Nm; ++j)
      {
        SymInvertChol(mpl[j]);
        SymInvertChol(mpz[j]);
      }
    }

    tmp = dtime() - t0;
    std::cout << "Matriplex invert time = " << tmp << " s\n";
    std::cout << "SMatrix / Matriplex = " << tsm/tmp << "\n\n";

    // ================================================================

    if (COMPARE_INVERSES)
    {
      for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
          printf("%8f ", mul[0](i,j)); printf("\n");
      } printf("\n");

      for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
          printf("%8f ", mpl[0].At(i,j,0)); printf("\n");
      } printf("\n");
    }
  }

  return 0;
}
