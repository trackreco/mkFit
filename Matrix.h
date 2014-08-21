#ifndef _matrix_
#define _matrix_

#include "Math/SMatrix.h"

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

#include "MatriplexSymNT.h"

const idx_t M = 6;

typedef Matriplex<float, M, M>   MPlexMM;
typedef Matriplex<float, M, 1>   MPlexMV;
typedef MatriplexSym<float, M>   MPlexSS;

//==============================================================================

#include <random>

extern std::default_random_engine            g_gen;
extern std::normal_distribution<float>       g_gaus;
extern std::uniform_real_distribution<float> g_unif;

#ifndef NO_ROOT

typedef double Double_t;

namespace TMath
{
   inline Double_t Pi()       { return 3.14159265358979323846; }
   inline Double_t TwoPi()    { return 2.0 * Pi(); }
   inline Double_t PiOver2()  { return Pi() / 2.0; }
   inline Double_t PiOver4()  { return Pi() / 4.0; }
   inline Double_t InvPi()    { return 1.0 / Pi(); }
   inline Double_t RadToDeg() { return 180.0 / Pi(); }
   inline Double_t DegToRad() { return Pi() / 180.0; }
   inline Double_t Sqrt2()    { return 1.4142135623730950488016887242097; }
}

#else

#include "TMath.h"

#endif

#endif
