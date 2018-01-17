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

typedef ROOT::Math::SMatrix<float,2> SMatrix22;
typedef ROOT::Math::SMatrix<float,2,2,ROOT::Math::MatRepSym<float,2> >    SMatrixSym22;
typedef ROOT::Math::SVector<float,2> SVector2;

typedef ROOT::Math::SMatrix<float,3,6> SMatrix36;
typedef ROOT::Math::SMatrix<float,6,3> SMatrix63;

typedef ROOT::Math::SMatrix<float,2,6> SMatrix26;
typedef ROOT::Math::SMatrix<float,6,2> SMatrix62;

// should work with any SMatrix
template<typename Matrix>
void dumpMatrix(Matrix m)
{
  for (int r=0;r<m.kRows;++r) {
    for (int c=0;c<m.kCols;++c) {
      std::cout << std::setw(12) << m.At(r,c) << " ";
    }
    std::cout << std::endl;
  }
}

template <typename Matrix>
inline void diagonalOnly(Matrix& m)
{
  for (int r=0; r<m.kRows; r++) {
    for (int c=0; c<m.kCols; c++) {
      if (r!=c) m[r][c] = 0.f;
    }
  }
}

//==============================================================================

// This should go elsewhere, eventually.

template<class T, class Compare> inline
constexpr const T clamp( const T v, const T lo, const T hi, Compare comp )
{
  return comp(v, lo) ? lo : comp(hi, v) ? hi : v;
}

template<class T> inline
constexpr const T clamp( const T v, const T lo, const T hi )
{
  return clamp( v, lo, hi, std::less<T>() );
}

#include <sys/time.h>

inline double dtime()
{
    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday(&mytime,(struct timezone*)0);
    tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
    return( tseconds );
}

CUDA_CALLABLE
inline float hipo(float x, float y)
{
  return std::sqrt(x*x + y*y);
}

CUDA_CALLABLE
inline void sincos4(const float x, float& sin, float& cos)
{
   // Had this writen with explicit division by factorial.
   // The *whole* fitting test ran like 2.5% slower on MIC, sigh.

   const float x2 = x*x;
   cos  = 1.f - 0.5f*x2 + 0.04166667f*x2*x2;
   sin  = x - 0.16666667f*x*x2;
}

//==============================================================================

// This ifdef needs to be changed to something like "use matriplex" and/or
// "is icc" as we can only do vectorization with icc now.

#ifdef USE_MATRIPLEX

  #ifdef __INTEL_COMPILER
    #define ASSUME_ALIGNED(a, b) __assume_aligned(a, b)
  #else
    #define ASSUME_ALIGNED(a, b) a = static_cast<decltype(a)>(__builtin_assume_aligned(a, b))
  #endif

  #include "Matriplex/MatriplexSym.h"

  constexpr Matriplex::idx_t NN =  MPT_SIZE; // "Length" of MPlex.

  constexpr Matriplex::idx_t LL =  6; // Dimension of large/long  MPlex entities
  constexpr Matriplex::idx_t HH =  3; // Dimension of small/short MPlex entities

  typedef Matriplex::Matriplex<float, LL, LL, NN>   MPlexLL;
  typedef Matriplex::Matriplex<float, LL,  1, NN>   MPlexLV;
  typedef Matriplex::MatriplexSym<float, LL,  NN>   MPlexLS;

  typedef Matriplex::Matriplex<float, HH, HH, NN>   MPlexHH;
  typedef Matriplex::Matriplex<float, HH,  1, NN>   MPlexHV;
  typedef Matriplex::MatriplexSym<float, HH,  NN>   MPlexHS;

  typedef Matriplex::Matriplex<float, 2,  2, NN>    MPlex22;
  typedef Matriplex::Matriplex<float, 2,  1, NN>    MPlex2V;
  typedef Matriplex::MatriplexSym<float,  2, NN>    MPlex2S;

  typedef Matriplex::Matriplex<float, LL, HH, NN>   MPlexLH;
  typedef Matriplex::Matriplex<float, HH, LL, NN>   MPlexHL;

  typedef Matriplex::Matriplex<float, LL,  2, NN>   MPlexL2;

  typedef Matriplex::Matriplex<float, 1, 1, NN>     MPlexQF;
  typedef Matriplex::Matriplex<int,   1, 1, NN>     MPlexQI;

  typedef Matriplex::Matriplex<bool,  1, 1, NN>     MPlexQB;

#endif

//==============================================================================

#include <random>

extern std::default_random_engine            g_gen;
extern std::normal_distribution<float>       g_gaus;
extern std::uniform_real_distribution<float> g_unif;

#endif
