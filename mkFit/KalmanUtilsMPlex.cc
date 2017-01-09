#include "KalmanUtilsMPlex.h"
#include "PropagationMPlex.h"

//#define DEBUG
#include "Debug.h"

#include "KalmanUtilsMPlex.icc"

namespace
{
  using idx_t = Matriplex::idx_t;

inline
void MultResidualsAdd(const MPlexLH& A,
		      const MPlexLV& B,
		      const MPlex2V& C,
		            MPlexLV& D)
{
   // outPar = psPar + kalmanGain*(dPar)
   //   D    =   B         A         C
   // where right half of kalman gain is 0 

   // XXX Regenerate with a script.

   MultResidualsAdd_imp(A, B, C, D, 0, NN);
}

inline
void MultResidualsAdd(const MPlexL2& A,
		      const MPlexLV& B,
		      const MPlex2V& C,
		            MPlexLV& D)
{
   // outPar = psPar + kalmanGain*(dPar)
   //   D    =   B         A         C
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
      // generate loop (can also write it manually this time, it's not much)
      d[0 * N + n] = b[0 * N + n] + a[ 0 * N + n] * c[0 * N + n] + a[ 1 * N + n] * c[1 * N + n];
      d[1 * N + n] = b[1 * N + n] + a[ 2 * N + n] * c[0 * N + n] + a[ 3 * N + n] * c[1 * N + n];
      d[2 * N + n] = b[2 * N + n] + a[ 4 * N + n] * c[0 * N + n] + a[ 5 * N + n] * c[1 * N + n];
      d[3 * N + n] = b[3 * N + n] + a[ 6 * N + n] * c[0 * N + n] + a[ 7 * N + n] * c[1 * N + n];
      d[4 * N + n] = b[4 * N + n] + a[ 8 * N + n] * c[0 * N + n] + a[ 9 * N + n] * c[1 * N + n];
      d[5 * N + n] = b[5 * N + n] + a[10 * N + n] * c[0 * N + n] + a[11 * N + n] * c[1 * N + n];
   }
}

//------------------------------------------------------------------------------

inline
void Chi2Similarity(const MPlex2V& A,//resPar
		    const MPlex2S& C,//resErr
                          MPlexQF& D)//outChi2
{

   // outChi2 = (resPar) * resErr * (resPar)
   //   D     =    A      *    C   *      A

   typedef float T;
   const idx_t N = NN;

   const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
   const T *c = C.fArray; ASSUME_ALIGNED(c, 64);
         T *d = D.fArray; ASSUME_ALIGNED(d, 64);

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      // generate loop (can also write it manually this time, it's not much)
      d[0 * N + n] = c[0 * N + n]*a[0 * N + n]*a[0 * N + n]
                   + c[2 * N + n]*a[1 * N + n]*a[1 * N + n] 
               + 2*( c[1 * N + n]*a[1 * N + n]*a[0 * N + n]);
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

inline
void AddIntoUpperLeft2x2(const MPlexLS& A, const MPlexHS& B, MPlex2S& C)
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
   }
}

//------------------------------------------------------------------------------

inline
void SubtractFirst3(const MPlexHV& A, const MPlexLV& B, MPlexHV& C)
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
      c[0*N+n] = a[0*N+n] - b[0*N+n];
      c[1*N+n] = a[1*N+n] - b[1*N+n];
      c[2*N+n] = a[2*N+n] - b[2*N+n];
   }
}

inline
void SubtractFirst2(const MPlexHV& A, const MPlexLV& B, MPlex2V& C)
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
      c[0*N+n] = a[0*N+n] - b[0*N+n];
      c[1*N+n] = a[1*N+n] - b[1*N+n];
   }
}

//==============================================================================

inline
void ProjectResErr(const MPlexQF& A00,
		   const MPlexQF& A01,
		   const MPlexHS& B, 
		         MPlexHH& C)
{
  // C = A * B, C is 3x3, A is 3x3 , B is 3x3 sym

  // Based on script generation and adapted to custom sizes.

  typedef float T;
  const idx_t N = NN;
  
  const T *a00 = A00.fArray; ASSUME_ALIGNED(a00, 64);
  const T *a01 = A01.fArray; ASSUME_ALIGNED(a01, 64);
  const T *b   = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c   = C.fArray; ASSUME_ALIGNED(c, 64);

#pragma simd
   for (int n = 0; n < N; ++n)
   {
      c[ 0*N+n] = a00[n]*b[ 0*N+n] + a01[n]*b[ 1*N+n];
      c[ 1*N+n] = a00[n]*b[ 1*N+n] + a01[n]*b[ 2*N+n];
      c[ 2*N+n] = a00[n]*b[ 3*N+n] + a01[n]*b[ 4*N+n];
      c[ 3*N+n] = b[ 3*N+n];
      c[ 4*N+n] = b[ 4*N+n];
      c[ 5*N+n] = b[ 5*N+n];
      c[ 6*N+n] = a01[n]*b[ 0*N+n] - a00[n]*b[ 1*N+n];
      c[ 7*N+n] = a01[n]*b[ 1*N+n] - a00[n]*b[ 2*N+n];
      c[ 8*N+n] = a01[n]*b[ 3*N+n] - a00[n]*b[ 4*N+n];
   }
}

inline
void ProjectResErrTransp(const MPlexQF& A00,
			 const MPlexQF& A01,
			 const MPlexHH& B, 
			       MPlex2S& C)
{
  // C = A * B, C is 3x3 sym, A is 3x3 , B is 3x3

  // Based on script generation and adapted to custom sizes.
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a00 = A00.fArray; ASSUME_ALIGNED(a00, 64);
  const T *a01 = A01.fArray; ASSUME_ALIGNED(a01, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#pragma simd
   for (int n = 0; n < N; ++n)
   {
      c[ 0*N+n] = b[ 0*N+n]*a00[n] + b[ 1*N+n]*a01[n];
      c[ 1*N+n] = b[ 3*N+n]*a00[n] + b[ 4*N+n]*a01[n];
      c[ 2*N+n] = b[ 5*N+n];
   }
}

#ifndef CCSCOORD
inline
void CCSErr(const MPlexLL& A, const MPlexLS& B, MPlexLL& C)
{
  // C = A * B, C is 6x6, A is 6x6 , B is 6x6 sym
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "CCSErr.ah"
}

inline
void CCSErrTransp(const MPlexLL& A, const MPlexLL& B, MPlexLS& C)
{
  // C = A * B, C is sym, A is 6x6 , B is 6x6
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "CCSErrTransp.ah"
}

inline
void CartesianErr(const MPlexLL& A, const MPlexLS& B, MPlexLL& C)
{
  // C = A * B, C is 6x6, A is 6x6 , B is 6x6 sym
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "CartesianErr.ah"
}

inline
void CartesianErrTransp(const MPlexLL& A, const MPlexLL& B, MPlexLS& C)
{
  // C = A * B, C is sym, A is 6x6 , B is 6x6
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "CartesianErrTransp.ah"
}
#endif

inline
void RotateResidulsOnTangentPlane(const MPlexQF& R00,//r00
				  const MPlexQF& R01,//r01
				  const MPlexHV& A  ,//res_glo
				        MPlex2V& B  )//res_loc
{
  RotateResidulsOnTangentPlane_impl(R00, R01, A, B, 0, NN);
#if 0
   // res_loc = rotT * res_glo
   //   B     =  R   *    A   

   typedef float T;
   const idx_t N = NN;

   const T *a   = A.fArray;   ASSUME_ALIGNED(a, 64);
   const T *r00 = R00.fArray; ASSUME_ALIGNED(r00, 64);
   const T *r01 = R01.fArray; ASSUME_ALIGNED(r01, 64);
         T *b   = B.fArray;   ASSUME_ALIGNED(b, 64);

#pragma simd
   for (idx_t n = 0; n < N; ++n)
   {
      b[0 * N + n] =  r00[0 * N + n]*a[0 * N + n] + r01[0 * N + n]*a[1 * N + n];
      b[1 * N + n] =  a[2 * N + n];
   }
#endif
}

inline
void KalmanHTG(const MPlexQF& A00,
	       const MPlexQF& A01,
	       const MPlex2S& B  ,
	             MPlexHH& C  )
{

   // HTG  = rot * res_loc
   //   C  =  A  *    B   

   // Based on script generation and adapted to custom sizes.

   typedef float T;
   const idx_t N = NN;

   const T *a00 = A00.fArray; ASSUME_ALIGNED(a00, 64);
   const T *a01 = A01.fArray; ASSUME_ALIGNED(a01, 64);
   const T *b   = B.fArray;   ASSUME_ALIGNED(b, 64);
         T *c   = C.fArray;   ASSUME_ALIGNED(c, 64);

#pragma simd
   for (int n = 0; n < N; ++n)
   {
      c[ 0*N+n] = a00[n]*b[ 0*N+n];
      c[ 1*N+n] = a00[n]*b[ 1*N+n];
      c[ 2*N+n] = 0.;
      c[ 3*N+n] = a01[n]*b[ 0*N+n];
      c[ 4*N+n] = a01[n]*b[ 1*N+n];
      c[ 5*N+n] = 0.;
      c[ 6*N+n] = b[ 1*N+n];
      c[ 7*N+n] = b[ 2*N+n];
      c[ 8*N+n] = 0.;
   }
}

inline
void KalmanGain(const MPlexLS& A, const MPlexHH& B, MPlexLH& C)
{
  // C = A * B, C is 6x3, A is 6x6 sym , B is 3x3

  typedef float T;
  const idx_t N = NN;

  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#pragma simd
   for (int n = 0; n < N; ++n)
   {
      c[ 0*N+n] = a[ 0*N+n]*b[ 0*N+n] + a[ 1*N+n]*b[ 3*N+n] + a[ 3*N+n]*b[ 6*N+n];
      c[ 1*N+n] = a[ 0*N+n]*b[ 1*N+n] + a[ 1*N+n]*b[ 4*N+n] + a[ 3*N+n]*b[ 7*N+n];
      c[ 2*N+n] = 0;
      c[ 3*N+n] = a[ 1*N+n]*b[ 0*N+n] + a[ 2*N+n]*b[ 3*N+n] + a[ 4*N+n]*b[ 6*N+n];
      c[ 4*N+n] = a[ 1*N+n]*b[ 1*N+n] + a[ 2*N+n]*b[ 4*N+n] + a[ 4*N+n]*b[ 7*N+n];
      c[ 5*N+n] = 0;
      c[ 6*N+n] = a[ 3*N+n]*b[ 0*N+n] + a[ 4*N+n]*b[ 3*N+n] + a[ 5*N+n]*b[ 6*N+n];
      c[ 7*N+n] = a[ 3*N+n]*b[ 1*N+n] + a[ 4*N+n]*b[ 4*N+n] + a[ 5*N+n]*b[ 7*N+n];
      c[ 8*N+n] = 0;
      c[ 9*N+n] = a[ 6*N+n]*b[ 0*N+n] + a[ 7*N+n]*b[ 3*N+n] + a[ 8*N+n]*b[ 6*N+n];
      c[10*N+n] = a[ 6*N+n]*b[ 1*N+n] + a[ 7*N+n]*b[ 4*N+n] + a[ 8*N+n]*b[ 7*N+n];
      c[11*N+n] = 0;
      c[12*N+n] = a[10*N+n]*b[ 0*N+n] + a[11*N+n]*b[ 3*N+n] + a[12*N+n]*b[ 6*N+n];
      c[13*N+n] = a[10*N+n]*b[ 1*N+n] + a[11*N+n]*b[ 4*N+n] + a[12*N+n]*b[ 7*N+n];
      c[14*N+n] = 0;
      c[15*N+n] = a[15*N+n]*b[ 0*N+n] + a[16*N+n]*b[ 3*N+n] + a[17*N+n]*b[ 6*N+n];
      c[16*N+n] = a[15*N+n]*b[ 1*N+n] + a[16*N+n]*b[ 4*N+n] + a[17*N+n]*b[ 7*N+n];
      c[17*N+n] = 0;
   }
}

void KalmanGain(const MPlexLS& A, const MPlex2S& B, MPlexL2& C)
{
  // C = A * B, C is 6x2, A is 6x6 sym , B is 2x2

  typedef float T;
  const idx_t N = NN;

  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "KalmanGain62.ah"
}

inline
void KHMult(const MPlexLH& A, 
	    const MPlexQF& B00,
	    const MPlexQF& B01,
	          MPlexLL& C)
{
  // C = A * B, C is 6x6, A is 6x3 , B is 3x6
  KHMult_imp(A, B00, B01, C, 0, NN);
}


inline
void KHC(const MPlexLL& A, const MPlexLS& B, MPlexLS& C)
{
  // C = A * B, C is 6x6, A is 6x6 , B is 6x6 sym

  typedef float T;
  const idx_t N = NN;

  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "KHC.ah"
}

inline
void KHC(const MPlexL2& A, const MPlexLS& B, MPlexLS& C)
{
  // C = A * B, C is 6x6 sym, A is 6x2 , B is 6x6 sym
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "K62HC.ah"
}

#ifndef CCSCOORD
inline
void ConvertToCCS(const MPlexLV& A, MPlexLV& B, MPlexLL& C)
{
  ConvertToCCS_imp(A, B, C, 0, NN);
}

inline
void ConvertToCartesian(const MPlexLV& A, MPlexLV& B, MPlexLL& C)
{
  ConvertToCartesian_imp(A, B, C, 0, NN); 
}
#endif

// //Warning: MultFull is not vectorized, use only for testing!
// template<typename T1, typename T2, typename T3>
// void MultFull(const T1& A, int nia, int nja, const T2& B, int nib, int njb, T3& C, int nic, int njc)
// {
// #ifdef DEBUG
//   assert(nja==nib);
//   assert(nia==nic);
//   assert(njb==njc);
// #endif
//   for (int n = 0; n < NN; ++n)
//     {
//       for (int i = 0; i < nia; ++i) {
// 	for (int j = 0; j < njb; ++j) {
// 	  C(n,i,j) = 0.;
// 	  for (int k = 0; k < nja; ++k) C(n,i,j) += A.ConstAt(n,i,k)*B.ConstAt(n,k,j);
// 	}
//       }
//     }
// }

// //Warning: MultTranspFull is not vectorized, use only for testing!
// // (careful about which one is transposed, I think rows and cols are swapped and the one that is transposed is A)
// template<typename T1, typename T2, typename T3>
// void MultTranspFull(const T1& A, int nia, int nja, const T2& B, int nib, int njb, T3& C, int nic, int njc)
// {
// #ifdef DEBUG
//   assert(nja==njb);
//   assert(nia==nic);
//   assert(nib==njc);
// #endif
//   for (int n = 0; n < NN; ++n)
//     {
//       for (int i = 0; i < nia; ++i) {
// 	for (int j = 0; j < nib; ++j) {
// 	  C(n,i,j) = 0.;
// 	  for (int k = 0; k < nja; ++k) C(n,i,j) += A.ConstAt(n,i,k)*B.ConstAt(n,j,k);
// 	}
//       }
//     }
// }

}


//==============================================================================
// updateParametersMPlex
//==============================================================================

void updateParametersMPlex(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                                 MPlexLS &outErr,       MPlexLV& outPar,
                           const int      N_proc)
{
  // const idx_t N = psErr.N;
  // Assert N-s of all parameters are the same.

  // Temporaries -- this is expensive -- should have them allocated outside and reused.
  // Can be passed in in a struct, see above.

  // updateParametersContext ctx;
  //assert((long long)(&updateCtx.propErr.fArray[0]) % 64 == 0);

  MPlexLS propErr;
  MPlexLV propPar;
  // do a full propagation step to correct for residual distance from the hit radius - need the charge for this
  if (Config::useCMSGeom) {
    propagateHelixToRMPlex(psErr,  psPar, inChg,  msPar, propErr, propPar, N_proc);
  } else {
    propErr = psErr;
    propPar = psPar;
  }

#ifdef DEBUG
  {
    dmutex_guard;
    printf("propPar:\n");
    for (int i = 0; i < 6; ++i) { 
      printf("%8f ", propPar.ConstAt(0,0,i)); printf("\n");
    } printf("\n");
    printf("msPar:\n");
    for (int i = 0; i < 3; ++i) { 
      printf("%8f ", msPar.ConstAt(0,0,i)); printf("\n");
    } printf("\n");
    printf("propErr:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", propErr.At(0,i,j)); printf("\n");
    } printf("\n");
    printf("msErr:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", msErr.ConstAt(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  // Rotate global point on tangent plane to cylinder
  // Tangent point is half way between hit and propagate position

  // Rotation matrix
  //  rotT00  0  rotT01
  //  rotT01  0 -rotT00
  //     0    1    0
  // Minimize temporaries: only two float are needed!

  MPlexQF rotT00;
  MPlexQF rotT01;
#pragma simd
  for (int n = 0; n < NN; ++n) {
    float r = hipo(msPar.ConstAt(n, 0, 0), msPar.ConstAt(n, 1, 0));
    rotT00.At(n, 0, 0) = -(msPar.ConstAt(n, 1, 0)+propPar.ConstAt(n, 1, 0))/(2*r);
    rotT01.At(n, 0, 0) =  (msPar.ConstAt(n, 0, 0)+propPar.ConstAt(n, 0, 0))/(2*r);
  }

  MPlexHV res_glo;   //position residual in global coordinates
  SubtractFirst3(msPar, propPar, res_glo);
  
  MPlexHS resErr_glo;//covariance sum in global position coordinates
  AddIntoUpperLeft3x3(propErr, msErr, resErr_glo);

  MPlex2V res_loc;   //position residual in local coordinates
  RotateResidulsOnTangentPlane(rotT00,rotT01,res_glo,res_loc);
  MPlex2S resErr_loc;//covariance sum in local position coordinates
  MPlexHH tempHH;
  ProjectResErr      (rotT00, rotT01, resErr_glo, tempHH);
  ProjectResErrTransp(rotT00, rotT01, tempHH, resErr_loc);

#ifdef DEBUG
  {
    dmutex_guard;
    printf("resErr:\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr_loc.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  //invert the 2x2 matrix
  Matriplex::InvertCramerSym(resErr_loc);

#ifndef CCSCOORD
  // Move to CCS coordinates: (x,y,z,1/pT,phi,theta)

  MPlexLV propPar_ccs;// propagated parameters in CCS coordinates
  MPlexLL jac_ccs;    // jacobian from cartesian to CCS
  ConvertToCCS(propPar,propPar_ccs,jac_ccs);

  MPlexLL tempLL;
  CCSErr      (jac_ccs, propErr, tempLL);
  CCSErrTransp(jac_ccs, tempLL, propErr);// propErr is now propagated errors in CCS coordinates
#endif

  // Kalman update in CCS coordinates

  MPlexLH K;           // kalman gain, fixme should be L2
  KalmanHTG(rotT00, rotT01, resErr_loc, tempHH); // intermediate term to get kalman gain (H^T*G)
  KalmanGain(propErr, tempHH, K);

#ifdef CCSCOORD
  MultResidualsAdd(K, propPar, res_loc, outPar);// propPar_ccs is now the updated parameters in CCS coordinates
  MPlexLL tempLL;
#else
  MultResidualsAdd(K, propPar_ccs, res_loc, propPar_ccs);// propPar_ccs is now the updated parameters in CCS coordinates
#endif


  KHMult(K, rotT00, rotT01, tempLL);
  KHC(tempLL, propErr, outErr);
  outErr.Subtract(propErr, outErr);// outErr is in CCS coordinates now

#ifndef CCSCOORD
  // Go back to cartesian coordinates

  // jac_ccs is now the jacobian from CCS to cartesian
  ConvertToCartesian(propPar_ccs, outPar, jac_ccs);
  CartesianErr      (jac_ccs, outErr, tempLL);
  CartesianErrTransp(jac_ccs, tempLL, outErr);// outErr is in cartesian coordinates now
#endif

#ifdef DEBUG
  {
    dmutex_guard;
    printf("res_glo:\n");
    for (int i = 0; i < 3; ++i) {
        printf("%8f ", res_glo.At(0,i,0));
    } printf("\n");
    printf("res_loc:\n");
    for (int i = 0; i < 2; ++i) {
        printf("%8f ", res_loc.At(0,i,0));
    } printf("\n");
    printf("resErr_loc (Inv):\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr_loc.At(0,i,j)); printf("\n");
    } printf("\n");
#ifndef CCSCOORD
    printf("jac_ccs:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", jac_ccs.At(0,i,j)); printf("\n");
    } printf("\n");
#endif
    printf("K:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", K.At(0,i,j)); printf("\n");
    } printf("\n");
    printf("outPar:\n");
    for (int i = 0; i < 6; ++i) {
      printf("%8f  ", outPar.At(0,i,0));
    } printf("\n");
    printf("outErr:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", outErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif
}

void computeChi2MPlex(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
                      const MPlexHS &msErr,  const MPlexHV& msPar,
                            MPlexQF& outChi2,
                      const int      N_proc)
{

  // const idx_t N = psErr.N;
  // Assert N-s of all parameters are the same.

  // Temporaries -- this is expensive -- should have them allocated outside and reused.
  // Can be passed in in a struct, see above.

  // updateParametersContext ctx;
  //assert((long long)(&updateCtx.propErr.fArray[0]) % 64 == 0);

  MPlexLS propErr;
  MPlexLV propPar;
  // do a full propagation step to correct for residual distance from the hit radius - need the charge for this
  if (Config::useCMSGeom) {
    propagateHelixToRMPlex(psErr,  psPar, inChg,  msPar, propErr, propPar, N_proc);
  } else {
    propErr = psErr;
    propPar = psPar;
  }

#ifdef DEBUG
  {
    dmutex_guard;
    printf("propPar:\n");
    for (int i = 0; i < 6; ++i) { 
      printf("%8f ", propPar.ConstAt(0,0,i)); printf("\n");
    } printf("\n");
    printf("propErr:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", propErr.At(0,i,j)); printf("\n");
    } printf("\n");
    printf("msPar:\n");
    for (int i = 0; i < 3; ++i) {
      printf("%8f ", msPar.ConstAt(0,0,i)); printf("\n");
    } printf("\n");
    printf("msErr:\n");
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", msErr.ConstAt(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  // Rotate global point on tangent plane to cylinder
  // Tangent point is half way between hit and propagate position

  // Rotation matrix
  //  rotT00  0  rotT01
  //  rotT01  0 -rotT00
  //     0    1    0
  // Minimize temporaries: only two float are needed!

  MPlexQF rotT00;
  MPlexQF rotT01;
  for (int n = 0; n < NN; ++n) {
    const float r = hipo(msPar.ConstAt(n, 0, 0), msPar.ConstAt(n, 1, 0));
    rotT00.At(n, 0, 0) = -(msPar.ConstAt(n, 1, 0)+propPar.ConstAt(n, 1, 0))/(2*r);
    rotT01.At(n, 0, 0) =  (msPar.ConstAt(n, 0, 0)+propPar.ConstAt(n, 0, 0))/(2*r);
  }

  MPlexHV res_glo;   //position residual in global coordinates
  SubtractFirst3(msPar, propPar, res_glo);
  
  MPlexHS resErr_glo;//covariance sum in global position coordinates
  AddIntoUpperLeft3x3(propErr, msErr, resErr_glo);

  MPlex2V res_loc;   //position residual in local coordinates
  RotateResidulsOnTangentPlane(rotT00,rotT01,res_glo,res_loc);
  MPlex2S resErr_loc;//covariance sum in local position coordinates
  MPlexHH tempHH;
  ProjectResErr      (rotT00, rotT01, resErr_glo, tempHH);
  ProjectResErrTransp(rotT00, rotT01, tempHH, resErr_loc);

#ifdef DEBUG
  {
    dmutex_guard;
    printf("resErr_loc:\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr_loc.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  //invert the 2x2 matrix
  Matriplex::InvertCramerSym(resErr_loc);

  //compute chi2
  Chi2Similarity(res_loc, resErr_loc, outChi2);

#ifdef DEBUG
  {
    dmutex_guard;
    printf("resErr_loc (Inv):\n");
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr_loc.At(0,i,j)); printf("\n");
    } printf("\n");
    printf("chi2: %8f\n", outChi2.At(0,0,0));
  }
#endif

}





void updateParametersEndcapMPlex(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
				 const MPlexHS &msErr,  const MPlexHV& msPar,
                                       MPlexLS &outErr,       MPlexLV& outPar,
				 const int      N_proc)
{
  // const idx_t N = psErr.N;
  // Assert N-s of all parameters are the same.

  // Temporaries -- this is expensive -- should have them allocated outside and reused.
  // Can be passed in in a struct, see above.

  // updateParametersContext ctx;
  //assert((long long)(&updateCtx.propErr.fArray[0]) % 64 == 0);

  MPlexLS propErr;
  MPlexLV propPar;
  // do a full propagation step to correct for residual distance from the hit radius - need the charge for this
  if (Config::useCMSGeom) {
    propagateHelixToZMPlex(psErr,  psPar, inChg,  msPar, propErr, propPar, N_proc);
  } else {
    propErr = psErr;
    propPar = psPar;
  }

#ifdef DEBUG
  {
    dmutex_guard;
  printf("updateParametersEndcapMPlex\n");
    printf("propPar:\n");
    for (int i = 0; i < 6; ++i) {
      printf("%8f ", propPar.ConstAt(0,0,i)); printf("\n");
    } printf("\n");
    printf("msPar:\n");
    for (int i = 0; i < 3; ++i) {
      printf("%8f ", msPar.ConstAt(0,0,i)); printf("\n");
    } printf("\n");
    printf("propErr:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", propErr.At(0,i,j)); printf("\n");
    } printf("\n");
    printf("msErr:\n");
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", msErr.ConstAt(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlex2V res;
  SubtractFirst2(msPar, propPar, res);

  MPlex2S resErr;
  AddIntoUpperLeft2x2(propErr, msErr, resErr);

#ifdef DEBUG
  {
    dmutex_guard;
    printf("resErr:\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  //invert the 2x2 matrix
  Matriplex::InvertCramerSym(resErr);

  MPlexL2 K;
  KalmanGain(propErr, resErr, K);

  MultResidualsAdd(K, propPar, res, outPar);

  KHC(K, propErr, outErr);

#ifdef DEBUG
  {
    printf("outErr before subtract:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
	printf("%8f ", outErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  outErr.Subtract(propErr, outErr);

#ifdef DEBUG
  {
    dmutex_guard;
    printf("res:\n");
    for (int i = 0; i < 2; ++i) {
        printf("%8f ", res.At(0,i,0));
    } printf("\n");
    printf("resErr (Inv):\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr.At(0,i,j)); printf("\n");
    } printf("\n");
    printf("K:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", K.At(0,i,j)); printf("\n");
    } printf("\n");
    printf("outPar:\n");
    for (int i = 0; i < 6; ++i) {
      printf("%8f  ", outPar.At(0,i,0));
    } printf("\n");
    printf("outErr:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", outErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif
}

void computeChi2EndcapMPlex(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
			    const MPlexHS &msErr,  const MPlexHV& msPar,
                                  MPlexQF& outChi2,
			    const int      N_proc)
{
  // const idx_t N = psErr.N;
  // Assert N-s of all parameters are the same.

  // Temporaries -- this is expensive -- should have them allocated outside and reused.
  // Can be passed in in a struct, see above.

  // updateParametersContext ctx;
  //assert((long long)(&updateCtx.propErr.fArray[0]) % 64 == 0);

  MPlexLS propErr;
  MPlexLV propPar;
  // do a full propagation step to correct for residual distance from the hit radius - need the charge for this
  if (Config::useCMSGeom) {
    propagateHelixToZMPlex(psErr,  psPar, inChg,  msPar, propErr, propPar, N_proc);
  } else {
    propErr = psErr;
    propPar = psPar;
  }

#ifdef DEBUG
  {
    dmutex_guard;
    printf("computeChi2EndcapMPlex\n");
    printf("propPar:\n");
    for (int i = 0; i < 6; ++i) {
      printf("%8f ", propPar.ConstAt(0,0,i)); printf("\n");
    } printf("\n");
    printf("msPar:\n");
    for (int i = 0; i < 3; ++i) {
      printf("%8f ", msPar.ConstAt(0,0,i)); printf("\n");
    } printf("\n");
    printf("propErr:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", propErr.At(0,i,j)); printf("\n");
    } printf("\n");
    printf("msErr:\n");
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", msErr.ConstAt(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlex2V res;
  SubtractFirst2(msPar, propPar, res);

  MPlex2S resErr;
  AddIntoUpperLeft2x2(propErr, msErr, resErr);

#ifdef DEBUG
  {
    dmutex_guard;
    printf("res:\n");
    for (int i = 0; i < 2; ++i) {
        printf("%8f ", res.At(0,0,i)); printf("\n");
    } printf("\n");
    printf("resErr:\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  //invert the 2x2 matrix
  Matriplex::InvertCramerSym(resErr);

  //compute chi2
  Chi2Similarity(res, resErr, outChi2);

#ifdef DEBUG
  {
    dmutex_guard;
    printf("resErr_loc (Inv):\n");
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr.At(0,i,j)); printf("\n");
    } printf("\n");
    printf("chi2: %8f\n", outChi2.At(0,0,0));
  }
#endif

}
