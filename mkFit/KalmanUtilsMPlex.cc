#include "KalmanUtilsMPlex.h"
#include "PropagationMPlex.h"

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
      d[1 * N + n] = b[1 * N + n] + a[ 3 * N + n] * c[0 * N + n] + a[ 4 * N + n] * c[1 * N + n];
      d[2 * N + n] = b[2 * N + n] + a[ 6 * N + n] * c[0 * N + n] + a[ 7 * N + n] * c[1 * N + n];
      d[3 * N + n] = b[3 * N + n] + a[ 9 * N + n] * c[0 * N + n] + a[10 * N + n] * c[1 * N + n];
      d[4 * N + n] = b[4 * N + n] + a[12 * N + n] * c[0 * N + n] + a[13 * N + n] * c[1 * N + n];
      d[5 * N + n] = b[5 * N + n] + a[15 * N + n] * c[0 * N + n] + a[16 * N + n] * c[1 * N + n];
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
      d[0 * N + n] =    c[0 * N + n]*a[0 * N + n]*a[0 * N + n] + c[2 * N + n]*a[1 * N + n]*a[1 * N + n] + 2*( c[1 * N + n]*a[1 * N + n]*a[0 * N + n]);
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

inline
void PolarErr(const MPlexLL& A, const MPlexLS& B, MPlexLL& C)
{
  // C = A * B, C is 6x6, A is 6x6 , B is 6x6 sym
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "PolarErr.ah"
}

inline
void PolarErrTransp(const MPlexLL& A, const MPlexLL& B, MPlexLS& C)
{
  // C = A * B, C is sym, A is 6x6 , B is 6x6
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "PolarErrTransp.ah"
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

inline
void RotateResidulsOnTangentPlane(const MPlexQF& R00,//r00
				  const MPlexQF& R01,//r01
				  const MPlexHV& A  ,//res_glo
				        MPlex2V& B  )//res_loc
{

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
  // C = A * B, C is 6x3, A is 6x6 sym , B is 6x3
 
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

inline
void KHMult(const MPlexLH& A, 
	    const MPlexQF& B00,
	    const MPlexQF& B01,
	          MPlexLL& C)
{
  // C = A * B, C is 6x6, A is 6x3 , B is 3x6
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b00 = B00.fArray; ASSUME_ALIGNED(b00, 64);
  const T *b01 = B01.fArray; ASSUME_ALIGNED(b01, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#pragma simd
   for (int n = 0; n < N; ++n)
   {
      c[ 0*N+n] = a[ 0*N+n]*b00[n];
      c[ 1*N+n] = a[ 0*N+n]*b01[n];
      c[ 2*N+n] = a[ 1*N+n];
      c[ 3*N+n] = 0;
      c[ 4*N+n] = 0;
      c[ 5*N+n] = 0;
      c[ 6*N+n] = a[ 3*N+n]*b00[n];
      c[ 7*N+n] = a[ 3*N+n]*b01[n];
      c[ 8*N+n] = a[ 4*N+n];
      c[ 9*N+n] = 0;
      c[10*N+n] = 0;
      c[11*N+n] = 0;
      c[12*N+n] = a[ 6*N+n]*b00[n];
      c[13*N+n] = a[ 6*N+n]*b01[n];
      c[14*N+n] = a[ 7*N+n];
      c[15*N+n] = 0;
      c[16*N+n] = 0;
      c[17*N+n] = 0;
      c[18*N+n] = a[ 9*N+n]*b00[n];
      c[19*N+n] = a[ 9*N+n]*b01[n];
      c[20*N+n] = a[10*N+n];
      c[21*N+n] = 0;
      c[22*N+n] = 0;
      c[23*N+n] = 0;
      c[24*N+n] = a[12*N+n]*b00[n];
      c[25*N+n] = a[12*N+n]*b01[n];
      c[26*N+n] = a[13*N+n];
      c[27*N+n] = 0;
      c[28*N+n] = 0;
      c[29*N+n] = 0;
      c[30*N+n] = a[15*N+n]*b00[n];
      c[31*N+n] = a[15*N+n]*b01[n];
      c[32*N+n] = a[16*N+n];
      c[33*N+n] = 0;
      c[34*N+n] = 0;
      c[35*N+n] = 0;
   }
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
void ConvertToPolar(const MPlexLV& A, MPlexLV& B, MPlexLL& C)
{
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
        T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#pragma simd
  for (int n = 0; n < N; ++n)
  {
    float pt = getHypot(a[ 3*N+n], a[ 4*N+n]);
    float p2 = pt*pt + a[ 5*N+n]*a[ 5*N+n];
    //
    b[ 0*N+n] = a[ 0*N+n];
    b[ 1*N+n] = a[ 1*N+n];
    b[ 2*N+n] = a[ 2*N+n];
    b[ 3*N+n] = 1./pt;
    b[ 4*N+n] = getPhi(a[ 3*N+n], a[ 4*N+n]); //fixme: use trig approx
    b[ 5*N+n] = getTheta(pt, a[ 5*N+n]);
    //
    c[ 0*N+n] = 1.;
    c[ 1*N+n] = 0.;
    c[ 2*N+n] = 0.;
    c[ 3*N+n] = 0.;
    c[ 4*N+n] = 0.;
    c[ 5*N+n] = 0.;
    c[ 6*N+n] = 0.;
    c[ 7*N+n] = 1.;
    c[ 8*N+n] = 0.;
    c[ 9*N+n] = 0.;
    c[10*N+n] = 0.;
    c[11*N+n] = 0.;
    c[12*N+n] = 0.;
    c[13*N+n] = 0.;
    c[14*N+n] = 1.;
    c[15*N+n] = 0.;
    c[16*N+n] = 0.;
    c[17*N+n] = 0.;
    c[18*N+n] = 0.;
    c[19*N+n] = 0.;
    c[20*N+n] = 0.;
    c[21*N+n] = -a[ 3*N+n]/(pt*pt*pt);
    c[22*N+n] = -a[ 4*N+n]/(pt*pt*pt);
    c[23*N+n] = 0.;
    c[24*N+n] = 0.;
    c[25*N+n] = 0.;
    c[26*N+n] = 0.;
    c[27*N+n] = -a[ 4*N+n]/(pt*pt);
    c[28*N+n] =  a[ 3*N+n]/(pt*pt);
    c[29*N+n] = 0.;
    c[30*N+n] = 0.;
    c[31*N+n] = 0.;
    c[32*N+n] = 0.;
    c[33*N+n] =  a[ 3*N+n]*a[ 5*N+n]/(pt*p2);
    c[34*N+n] =  a[ 4*N+n]*a[ 5*N+n]/(pt*p2);
    c[35*N+n] = -pt/p2;
  }
}

inline
void ConvertToCartesian(const MPlexLV& A, MPlexLV& B, MPlexLL& C)
{
 
  typedef float T;
  const idx_t N = NN;
  
  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
        T *b = B.fArray; ASSUME_ALIGNED(b, 64);
        T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#pragma simd
  for (int n = 0; n < N; ++n)
  {
    float cosP = cos(a[ 4*N+n]); //fixme: use trig approx
    float sinP = sin(a[ 4*N+n]);
    float cosT = cos(a[ 5*N+n]);
    float sinT = sin(a[ 5*N+n]);
    //
    b[ 0*N+n] = a[ 0*N+n];
    b[ 1*N+n] = a[ 1*N+n];
    b[ 2*N+n] = a[ 2*N+n];
    b[ 3*N+n] = cosP/a[ 3*N+n];
    b[ 4*N+n] = sinP/a[ 3*N+n];
    b[ 5*N+n] = cosT/(sinT*a[ 3*N+n]);
    //
    c[ 0*N+n] = 1.;
    c[ 1*N+n] = 0.;
    c[ 2*N+n] = 0.;
    c[ 3*N+n] = 0.;
    c[ 4*N+n] = 0.;
    c[ 5*N+n] = 0.;
    c[ 6*N+n] = 0.;
    c[ 7*N+n] = 1.;
    c[ 8*N+n] = 0.;
    c[ 9*N+n] = 0.;
    c[10*N+n] = 0.;
    c[11*N+n] = 0.;
    c[12*N+n] = 0.;
    c[13*N+n] = 0.;
    c[14*N+n] = 1.;
    c[15*N+n] = 0.;
    c[16*N+n] = 0.;
    c[17*N+n] = 0.;
    c[18*N+n] = 0.;
    c[19*N+n] = 0.;
    c[20*N+n] = 0.;
    c[21*N+n] = -cosP/(a[ 3*N+n]*a[ 3*N+n]);
    c[22*N+n] = -sinP/a[ 3*N+n];
    c[23*N+n] = 0.;
    c[24*N+n] = 0.;
    c[25*N+n] = 0.;
    c[26*N+n] = 0.;
    c[27*N+n] = -sinP/(a[ 3*N+n]*a[ 3*N+n]);
    c[28*N+n] =  cosP/a[ 3*N+n];
    c[29*N+n] = 0.;
    c[30*N+n] = 0.;
    c[31*N+n] = 0.;
    c[32*N+n] = 0.;
    c[33*N+n] = -cosT/(sinT*a[ 3*N+n]*a[ 3*N+n]);
    c[34*N+n] = 0.;
    c[35*N+n] = -1./(sinT*sinT*a[ 3*N+n]);
  }
}


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

//#define DEBUG

void updateParametersMPlex(const MPlexLS &psErr,  const MPlexLV& psPar, const MPlexQI &inChg,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                                MPlexLS &outErr,       MPlexLV& outPar)
{
  // const idx_t N = psErr.N;
  // Assert N-s of all parameters are the same.

  // Temporaries -- this is expensive -- should have them allocated outside and reused.
  // Can be passed in in a struct, see above.

#ifdef DEBUG
  const bool dump = g_dump;
#endif

  // updateParametersContext ctx;
  //assert((long long)(&updateCtx.propErr.fArray[0]) % 64 == 0);

  MPlexLS propErr;
  MPlexLV propPar;
  // do a full propagation step to correct for residual distance from the hit radius - need the charge for this
  if (Config::useCMSGeom) {
    propagateHelixToRMPlex(psErr,  psPar, inChg,  msPar, propErr, propPar);
  } else {
    propErr = psErr;
    propPar = psPar;
  }

#ifdef DEBUG
  if (dump) {
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
  if (dump) {
    printf("resErr:\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr_loc.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  //invert the 2x2 matrix
  Matriplex::InvertCramerSym(resErr_loc);

  // Move to "polar" coordinates: (x,y,z,1/pT,phi,theta) [can we find a better name?]

  MPlexLV propPar_pol;// propagated parameters in "polar" coordinates
  MPlexLL jac_pol;    // jacobian from cartesian to "polar"
  ConvertToPolar(propPar,propPar_pol,jac_pol);

  MPlexLL tempLL;
  PolarErr      (jac_pol, propErr, tempLL);
  PolarErrTransp(jac_pol, tempLL, propErr);// propErr is now propagated errors in "polar" coordinates

  // Kalman update in "polar" coordinates

  MPlexLH K;           // kalman gain, fixme should be L2
  KalmanHTG(rotT00, rotT01, resErr_loc, tempHH); // intermediate term to get kalman gain (H^T*G)
  KalmanGain(propErr, tempHH, K);

  MultResidualsAdd(K, propPar_pol, res_loc, propPar_pol);// propPar_pol is not the updated parameters in "polar" coordinates

  KHMult(K, rotT00, rotT01, tempLL);
  KHC(tempLL, propErr, outErr);
  outErr.Subtract(propErr, outErr);// outErr is in "polar" coordinates now

  // Go back to cartesian coordinates

  // jac_pol is now the jacobian from "polar" to cartesian
  ConvertToCartesian(propPar_pol, outPar, jac_pol);
  CartesianErr      (jac_pol, outErr, tempLL);
  CartesianErrTransp(jac_pol, tempLL, outErr);// outErr is in cartesian coordinates now

#ifdef DEBUG
  if (dump) {
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
    printf("jac_pol:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", jac_pol.At(0,i,j)); printf("\n");
    } printf("\n");
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
                            MPlexQF& outChi2)
{

  // const idx_t N = psErr.N;
  // Assert N-s of all parameters are the same.

  // Temporaries -- this is expensive -- should have them allocated outside and reused.
  // Can be passed in in a struct, see above.

#ifdef DEBUG
  const bool dump = g_dump;
#endif

  // updateParametersContext ctx;
  //assert((long long)(&updateCtx.propErr.fArray[0]) % 64 == 0);

  MPlexLS propErr;
  MPlexLV propPar;
  // do a full propagation step to correct for residual distance from the hit radius - need the charge for this
  if (Config::useCMSGeom) {
    propagateHelixToRMPlex(psErr,  psPar, inChg,  msPar, propErr, propPar);
  } else {
    propErr = psErr;
    propPar = psPar;
  }

#ifdef DEBUG
  if (dump) {
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
  if (dump) {
    printf("resErr_loc:\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr_loc.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  //invert the 2x2 matrix
  Matriplex::InvertCramerSym(resErr_loc);

#ifdef DEBUG
  if (dump) {
    printf("resErr_loc (Inv):\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr_loc.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  //compute chi2
  Chi2Similarity(res_loc, resErr_loc, outChi2);
}
