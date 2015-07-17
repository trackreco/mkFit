#include "KalmanUtilsMPlex.h"

namespace
{
  using idx_t = Matriplex::idx_t;

inline
void MultResidualsAdd(const MPlexLH& A,
                      const MPlexLV& B,
                      const MPlexHV& C,
                            MPlexLV& D)
{
   // outPar = psPar + kalmanGain*(msPar-psPar)
   //   D    =   B         A         C  -  B
   // where right half of kalman gain is 0 

   // XXX Regenerate with a script.

   typedef float T;
   const idx_t N = NN;

   const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
   const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
   const T *c = C.fArray; ASSUME_ALIGNED(c, 64);
         T *d = D.fArray; ASSUME_ALIGNED(d, 64);

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      // manually subrtact into local vars -- 3 of them
      float x0 = c[0 * N + n] - b[0 * N + n];
      float x1 = c[1 * N + n] - b[1 * N + n];
      float x2 = c[2 * N + n] - b[2 * N + n];

      // generate loop (can also write it manually this time, it's not much)
      d[0 * N + n] = b[0 * N + n] + a[ 0 * N + n] * x0 + a[ 1 * N + n] * x1 + a[ 2 * N + n] * x2;
      d[1 * N + n] = b[1 * N + n] + a[ 3 * N + n] * x0 + a[ 4 * N + n] * x1 + a[ 5 * N + n] * x2;
      d[2 * N + n] = b[2 * N + n] + a[ 6 * N + n] * x0 + a[ 7 * N + n] * x1 + a[ 8 * N + n] * x2;
      d[3 * N + n] = b[3 * N + n] + a[ 9 * N + n] * x0 + a[10 * N + n] * x1 + a[11 * N + n] * x2;
      d[4 * N + n] = b[4 * N + n] + a[12 * N + n] * x0 + a[13 * N + n] * x1 + a[14 * N + n] * x2;
      d[5 * N + n] = b[5 * N + n] + a[15 * N + n] * x0 + a[16 * N + n] * x1 + a[17 * N + n] * x2;
   }
}

//------------------------------------------------------------------------------

inline
void Chi2Similarity(const MPlexHV& A,//msPar
		    const MPlexLV& B,//psPar
		    const MPlexHS& C,//resErr
		    MPlexQF& D)//outChi2
{

   // outChi2 = (msPar-psPar) * resErr * (msPar-psPar)
   //   D     =    A - B      *    C   *      A-B

   // XXX Regenerate with a script.

   typedef float T;
   const idx_t N = NN;

   const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
   const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
   const T *c = C.fArray; ASSUME_ALIGNED(c, 64);
         T *d = D.fArray; ASSUME_ALIGNED(d, 64);

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      // manually subrtact into local vars -- 3 of them
      float x0 = a[0 * N + n] - b[0 * N + n];
      float x1 = a[1 * N + n] - b[1 * N + n];
      float x2 = a[2 * N + n] - b[2 * N + n];

      /*
      std::cout << "x: " << x0 << ", " << x1 << ", " << x2 << std::endl;
      std::cout << "c0: " << c[0 * N + n] << ", " << c[1 * N + n] << ", " << c[3 * N + n] << std::endl;
      std::cout << "c0: " << c[1 * N + n] << ", " << c[2 * N + n] << ", " << c[4 * N + n] << std::endl;
      std::cout << "c0: " << c[3 * N + n] << ", " << c[4 * N + n] << ", " << c[5 * N + n] << std::endl;
      */

      // generate loop (can also write it manually this time, it's not much)
      d[0 * N + n] = c[0 * N + n]*x0*x0 + c[2 * N + n]*x1*x1 + c[5 * N + n]*x2*x2 + 2*( c[1 * N + n]*x1*x0 + c[3 * N + n]*x2*x0 + c[4 * N + n]*x1*x2  );

   }
}

//------------------------------------------------------------------------------

inline
void AddIntoUpperLeft3x3(const MPlexLS& A, const MPlexHS& B, MPlexHS& C)
{
   // The rest of matrix is left untouched.

   typedef float T;
   const idx_t N = NN;

   const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
   const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
         T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      c[0*N+n] = a[0*N+n] + b[0*N+n];
      c[1*N+n] = a[1*N+n] + b[1*N+n];
      c[2*N+n] = a[2*N+n] + b[2*N+n];
      c[3*N+n] = a[3*N+n] + b[3*N+n];
      c[4*N+n] = a[4*N+n] + b[4*N+n];
      c[5*N+n] = a[5*N+n] + b[5*N+n];
   }
}

//==============================================================================

void MultKalmanGain(const MPlexLS& A, const MPlexHS& B, MPlexLH& C)
{
  // C = A * B, C is 6x3, A is 6x6 sym, B is 3x3 sym

  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "upParam_MultKalmanGain.ah"

}

void simil_x_propErr(const MPlexHS& A, const MPlexLS& B, MPlexLL& C)
{
  // C = A * B, C is 6x6, A is 3x3 sym, B is 6x6 sym, yes we're cheating with the math by making a copy of 3x3 into 6x6 at the same time as doing the actual multiplication
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "upParam_simil_x_propErr.ah"
}

void propErrT_x_simil_propErr(const MPlexLH& A, const MPlexLS& B, MPlexLS& C)
{
  // C = A * B, C is 6x6 sym, A is 6x3, B is 6x6 sym
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "upParam_propErrT_x_simil_propErr.ah"
}

void kalmanGain_x_propErr(const MPlexLH& A, const MPlexLS& B, MPlexLS& C)
{
  // C = A * B, C is 6x6 sym, A is 6x6 , B is 6x6 sym
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "upParam_kalmanGain_x_propErr.ah"
}

}


//==============================================================================
// updateParametersMPlex
//==============================================================================

void updateParametersMPlex(const MPlexLS &psErr,  const MPlexLV& psPar,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                           MPlexLS &outErr,       MPlexLV& outPar)
{
  // const idx_t N = psErr.N;
  // Assert N-s of all parameters are the same.

  // Temporaries -- this is expensive -- should have them allocated outside and reused.
  // Can be passed in in a struct, see above.

  // Also: resErr could be 3x3, kalmanGain 6x3

#ifdef DEBUG
  const bool dump = g_dump;
#endif

  // updateParametersContext ctx;
  //assert((long long)(&updateCtx.propErr.fArray[0]) % 64 == 0);

  MPlexLS propErr;
  propErr = psErr;       // could use/overwrite psErr?
  propErr.AddNoiseIntoUpperLeft3x3(0.0); // e.g. ?

#ifdef DEBUG
  if (dump) {
    printf("propErr:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", propErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  // if (dump) {
  //   printf("msErr:\n");
  //   for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
  //       printf("%8f ", msErr.ConstAt(0,i,j)); printf("\n");
  //   } printf("\n");
  // }

  MPlexHS resErr;
  AddIntoUpperLeft3x3(propErr, msErr, resErr);
  // Do not really need to zero the rest ... it is not used.

#ifdef DEBUG
  if (dump) {
    printf("resErr:\n");
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", resErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  Matriplex::InvertCramerSym(resErr);
  // Matriplex::InvertCholeskySym(resErr);
  // resErr is now resErrInv
  // XXX Both could be done in one operation.

#ifdef DEBUG
  if (dump) {
    printf("resErrInv:\n");
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", resErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexLH kalmanGain;
  MultKalmanGain(propErr, resErr, kalmanGain);

#ifdef DEBUG
  if (dump) {
    printf("kalmanGain:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", kalmanGain.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  // outPar = psPar + kalmanGain*(msPar-psPar)[0-2] (last three are 0!)
  MultResidualsAdd(kalmanGain, psPar, msPar, outPar);

#ifdef DEBUG
  if (dump) {
    printf("outPar:\n");
    for (int i = 0; i < 6; ++i) {
      printf("%8f  ", outPar.At(0,i,0));
    } printf("\n");
  }
#endif

  // outErr = propErr - ROOT::Math::SimilarityT(propErr,resErrInv);
  //        = propErr - kalmanGain*propErr
  //
  // XXX Ideally would also subtract at the same time in auto generated code.

  MPlexLS outErrTemp;
  kalmanGain_x_propErr(kalmanGain, propErr, outErrTemp);
  outErr.Subtract(propErr, outErrTemp);

#ifdef DEBUG
  if (dump) {
    printf("outErr:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", outErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif
}


void computeChi2MPlex(const MPlexLS &psErr,  const MPlexLV& psPar,
		      const MPlexHS &msErr,  const MPlexHV& msPar,
		      MPlexQF& outChi2)
{

  // const idx_t N = psErr.N;
  // Assert N-s of all parameters are the same.

  // Temporaries -- this is expensive -- should have them allocated outside and reused.
  // Can be passed in in a struct, see above.

  // Also: resErr could be 3x3, kalmanGain 6x3

#ifdef DEBUG
  const bool dump = true;//g_dump;
#endif

  // updateParametersContext ctx;
  //assert((long long)(&updateCtx.propErr.fArray[0]) % 64 == 0);

  MPlexLS propErr;
  propErr = psErr;       // could use/overwrite psErr?
  propErr.AddNoiseIntoUpperLeft3x3(0.0); // e.g. ?

#ifdef DEBUG
  if (dump) {
    MPlexLV propState = psPar;
    printf("propState:\n");
    for (int i = 0; i < 6; ++i) { 
      printf("%8f ", propState.At(0,0,i));
    } printf("\n");
    printf("propErr:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", propErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  // if (dump) {
  //   printf("msErr:\n");
  //   for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
  //       printf("%8f ", msErr.ConstAt(0,i,j)); printf("\n");
  //   } printf("\n");
  // }

  MPlexHS resErr;
  AddIntoUpperLeft3x3(propErr, msErr, resErr);
  // Do not really need to zero the rest ... it is not used.

#ifdef DEBUG
  if (dump) {
    printf("resErr:\n");
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", resErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  Matriplex::InvertCramerSym(resErr);
  // Matriplex::InvertCholeskySym(resErr);
  // resErr is now resErrInv
  // XXX Both could be done in one operation.

#ifdef DEBUG
  if (dump) {
    printf("resErrInv:\n");
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", resErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  Chi2Similarity(msPar,psPar,resErr, outChi2);

}
