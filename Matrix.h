#ifndef _matrix_
#define _matrix_

#include "Math/SMatrix.h"
#include "Config.h"

typedef ROOT::Math::SMatrix<float,6,6,ROOT::Math::MatRepSym<float,6> >    SMatrixSym66;
typedef ROOT::Math::SMatrix<float,6> SMatrix66;
typedef ROOT::Math::SVector<float,6> SVector6;

typedef ROOT::Math::SMatrix<float,3> SMatrix33;
typedef ROOT::Math::SMatrix<float,3,3,ROOT::Math::MatRepSym<float,3> >    SMatrixSym33;
typedef ROOT::Math::SVector<float,3> SVector3;

typedef ROOT::Math::SMatrix<float,3,6> SMatrix36;
typedef ROOT::Math::SMatrix<float,6,3> SMatrix63;

// should work with any SMatrix
template<typename Matrix>
void dumpMatrix(Matrix m)
{
  for (unsigned int r=0;r<m.kRows;++r) {
    for (unsigned int c=0;c<m.kCols;++c) {
      std::cout << std::setw(12) << m.At(r,c) << " ";
    }
    std::cout << std::endl;
  }
}


//==============================================================================

// This should go elsewhere, eventually.

#include <sys/time.h>

inline double dtime()
{
    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday(&mytime,(struct timezone*)0);
    tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
    return( tseconds );
}

inline float hipo(float x, float y)
{
   return sqrt(x*x + y*y);
}

inline void sincos4(float x, float& sin, float& cos)
{
   // Had this writen with explicit division by factorial.
   // The *whole* fitting test ran like 2.5% slower on MIC, sigh.

   cos  = 1;
   sin  = x;   x *= x * 0.5f;
   cos -= x;   x *= x * 0.33333333f;
   sin -= x;   x *= x * 0.25f;
   cos += x;
}

//==============================================================================

// This ifdef needs to be changed to something like "use matriplex" and/or
// "is icc" as we can only do vectorization with icc now.

#ifdef USE_MATRIPLEX

  #ifdef __INTEL_COMPILER

    #define ASSUME_ALIGNED(a, b) __assume_aligned(a, b)

  #else

    template<typename T> inline void ASSUME_ALIGNED(T* a, int b) { a = (T*) __builtin_assume_aligned(a, b); }

  #endif

  #include "Matriplex/MatriplexSym.h"

  const Matriplex::idx_t NN =  MPT_SIZE; // "Length" of MPlex.

  const Matriplex::idx_t LL =  6; // Dimension of large/long  MPlex entities
  const Matriplex::idx_t HH =  3; // Dimension of small/short MPlex entities

  typedef Matriplex::Matriplex<float, LL, LL, NN>   MPlexLL;
  typedef Matriplex::Matriplex<float, LL,  1, NN>   MPlexLV;
  typedef Matriplex::MatriplexSym<float, LL,  NN>   MPlexLS;

  typedef Matriplex::Matriplex<float, HH, HH, NN>   MPlexHH;
  typedef Matriplex::Matriplex<float, HH,  1, NN>   MPlexHV;
  typedef Matriplex::MatriplexSym<float, HH,  NN>   MPlexHS;

  typedef Matriplex::Matriplex<float, LL, HH, NN>   MPlexLH;

  typedef Matriplex::Matriplex<float, 1, 1, NN>     MPlexQF;
  typedef Matriplex::Matriplex<int,   1, 1, NN>     MPlexQI;

#endif

//==============================================================================

#include <random>

extern std::default_random_engine            g_gen;
extern std::normal_distribution<float>       g_gaus;
extern std::uniform_real_distribution<float> g_unif;

// All debug printouts are ifdefed with DEBUG
// #define DEBUG

extern bool g_dump;

#endif
