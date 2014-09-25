#ifndef _matrix_
#define _matrix_

#ifndef MAX_HITS
#define MAX_HITS 10
#endif

#include "Math/SMatrix.h"

typedef ROOT::Math::SMatrix<float,6,6,ROOT::Math::MatRepSym<float,6> >    SMatrixSym66;
typedef ROOT::Math::SMatrix<float,6> SMatrix66;
typedef ROOT::Math::SVector<float,6> SVector6;

typedef ROOT::Math::SMatrix<float,3> SMatrix33;
typedef ROOT::Math::SMatrix<float,3,3,ROOT::Math::MatRepSym<float,3> >    SMatrixSym33;
typedef ROOT::Math::SVector<float,3> SVector3;

typedef ROOT::Math::SMatrix<float,3,6> SMatrix36;
typedef ROOT::Math::SMatrix<float,6,3> SMatrix63;

void dumpMatrix(SMatrix33& m);
void dumpMatrix(SMatrix36& m);
void dumpMatrix(SMatrix63& m);
void dumpMatrix(SMatrix66& m);
void dumpMatrix(SMatrixSym33& m);
void dumpMatrix(SMatrixSym66& m);

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

//==============================================================================

// This ifdef needs to be changed to something like "use matriplex" and/or
// "is icc" as we can only do vectorization with icc now.

#ifndef __APPLE__

#include "Matriplex/MatriplexSym.h"

#ifndef MPT_SIZE
#define MPT_SIZE 16
#endif

typedef int idx_t;

const idx_t NN =  MPT_SIZE; // "Length" of MPlex.

const idx_t LL =  6; // Dimension of large/long  MPlex entities
const idx_t HH =  3; // Dimension of small/short MPlex entities

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

#ifndef NUM_THREADS
#define NUM_THREADS 1
#endif

#ifndef THREAD_BINDING
#define THREAD_BINDING spread
#endif

//==============================================================================

#include <random>

extern std::default_random_engine            g_gen;
extern std::normal_distribution<float>       g_gaus;
extern std::uniform_real_distribution<float> g_unif;

// All debug printouts are ifdefed with DEBUG
// #define DEBUG

extern bool g_dump;

#ifdef NO_ROOT

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
