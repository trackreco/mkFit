#include "KalmanUtilsMPlex.h"
#include "PropagationMPlex.h"

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
      // manually subtract into local vars -- 3 of them
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
      d[0 * N + n] =    c[0 * N + n]*x0*x0 + c[2 * N + n]*x1*x1 + c[5 * N + n]*x2*x2 +
                    2*( c[1 * N + n]*x1*x0 + c[3 * N + n]*x2*x0 + c[4 * N + n]*x1*x2);
   }
}

inline
void Chi2Similarity(const MPlexHV& A,//resPar
		    const MPlex2S& C,//resErr
                          MPlexQF& D)//outChi2
{

   // outChi2 = (resPar) * resErr * (resPar)
   //   D     =    A      *    C   *      A

   // XXX Regenerate with a script.

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

template<typename T1, typename T2, typename T3>
void MultFull(const T1& A, int nia, int nja, const T2& B, int nib, int njb, T3& C, int nic, int njc)
{

  assert(nja==nib);
  assert(nia==nic);
  assert(njb==njc);

#pragma simd
  for (int n = 0; n < NN; ++n)
    {
      for (int i = 0; i < nia; ++i) {
	for (int j = 0; j < njb; ++j) {
	  C(n,i,j) = 0.;
	  for (int k = 0; k < nja; ++k) C(n,i,j) += A.ConstAt(n,i,k)*B.ConstAt(n,k,j);
	}
      }
    }
}

template<typename T1, typename T2, typename T3>
void MultTranspFull(const T1& A, int nia, int nja, const T2& B, int nib, int njb, T3& C, int nic, int njc)
{

  assert(nja==njb);
  assert(nia==nic);
  assert(nib==njc);

#pragma simd
  for (int n = 0; n < NN; ++n)
    {
      for (int i = 0; i < nia; ++i) {
	for (int j = 0; j < nib; ++j) {
	  C(n,i,j) = 0.;
	  for (int k = 0; k < nja; ++k) C(n,i,j) += A.ConstAt(n,i,k)*B.ConstAt(n,j,k);
	}
      }
    }
}

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

  // Also: resErr could be 3x3, kalmanGain 6x3

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

  MPlexHH rot;//we may not need this one
  MPlexHH rotT;
  MPlexHV res_glo;
  MPlexHS resErr_glo;
#pragma simd
  for (int n = 0; n < NN; ++n) {
    float r = hipo(msPar.ConstAt(n, 0, 0), msPar.ConstAt(n, 1, 0));
    rot.At(n, 0, 0) = -(msPar.ConstAt(n, 1, 0)+propPar.ConstAt(n, 1, 0))/(2*r);
    rot.At(n, 0, 1) = 0;
    rot.At(n, 0, 2) =  (msPar.ConstAt(n, 0, 0)+propPar.ConstAt(n, 0, 0))/(2*r);
    rot.At(n, 1, 0) = rot.ConstAt(n, 0, 2);
    rot.At(n, 1, 1) = 0;
    rot.At(n, 1, 2) = -rot.ConstAt(n, 0, 0);
    rot.At(n, 2, 0) = 0;
    rot.At(n, 2, 1) = 1;
    rot.At(n, 2, 2) = 0;
    //
    rotT.At(n, 0, 0) = rot.ConstAt(n, 0, 0);
    rotT.At(n, 0, 1) = rot.ConstAt(n, 1, 0);
    rotT.At(n, 0, 2) = rot.ConstAt(n, 2, 0);
    rotT.At(n, 1, 0) = rot.ConstAt(n, 0, 1);
    rotT.At(n, 1, 1) = rot.ConstAt(n, 1, 1);
    rotT.At(n, 1, 2) = rot.ConstAt(n, 2, 1);
    rotT.At(n, 2, 0) = rot.ConstAt(n, 0, 2);
    rotT.At(n, 2, 1) = rot.ConstAt(n, 1, 2);
    rotT.At(n, 2, 2) = rot.ConstAt(n, 2, 2);
    //
    for (int i = 0; i < 3; ++i) {
      res_glo.At(n, i, 0) = msPar.ConstAt(n, i, 0) - propPar.ConstAt(n, i, 0);
      for (int j = 0; j < 3; ++j) {
	resErr_glo.At(n, i, j) = msErr.ConstAt(n,i,j) + propErr.ConstAt(n,i,j);
      }
    }
  }

  MPlexHV res_loc;
  MultiplyGeneral(rotT,res_glo,res_loc);

#ifdef DEBUG
  if (dump) {
    printf("res_loc:\n");
    for (int i = 0; i < 3; ++i) {
        printf("%8f ", res_loc.At(0,i,0));
    } printf("\n");
  }
#endif

  MPlexHS resErr_loc;
  MPlexHH temp;
  MultFull      (rotT,3,3, resErr_glo,3,3, temp,3,3);
  MultTranspFull(rotT,3,3, temp,3,3, resErr_loc,3,3);

  MPlex2S resErr;
#pragma simd
  for (int n = 0; n < NN; ++n) {
    for (int i = 0; i < 2; ++i) {
      for (int j = i; j < 2; ++j) {
	resErr.At(n, i, j) = resErr_loc.ConstAt(n,i,j);
      }
    }
  }

#ifdef DEBUG
  if (dump) {
    printf("resErr:\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  Matriplex::InvertCramerSym(resErr);

#ifdef DEBUG
  if (dump) {
    printf("resErrInv:\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

#pragma simd
  for (int n = 0; n < NN; ++n) {
    for (int i = 0; i < 3; ++i) {
      for (int j = i; j < 3; ++j) {
	if (i==2||j==2) resErr_loc.At(n, i, j) = 0.;
	else resErr_loc.At(n, i, j) = resErr.ConstAt(n,i,j);
      }
    }
  }


#ifdef DEBUG
  if (dump) {
    printf("resErr_loc:\n");
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", resErr_loc.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexHH resErrTmp;
  MultFull(rot,3,3, resErr_loc,3,3, resErrTmp,3,3);

#ifdef DEBUG
  if (dump) {
    printf("resErrTmp:\n");
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", resErrTmp.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexLH resErrTmpLH;
#pragma simd
  for (int n = 0; n < NN; ++n) {
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 3; ++j) {
	if (i>2) resErrTmpLH.At(n, i, j) = 0.;
	else resErrTmpLH.At(n, i, j) = resErrTmp.ConstAt(n,i,j);
      }
    }
  }

#ifdef DEBUG
  if (dump) {
    printf("resErrTmpLH:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", resErrTmpLH.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexLV propPar_pol;
  MPlexLL jac_pol;
#pragma simd
  for (int n = 0; n < NN; ++n) {
    //propPar_pol
    for (int i = 0; i < 3; ++i) {
      propPar_pol.At(n, i, 0) = propPar.At(n, i, 0);
    }
    float pt = getHypot(propPar.At(n, 3, 0), propPar.At(n, 4, 0));
    propPar_pol.At(n, 3, 0) = 1./pt;
    propPar_pol.At(n, 4, 0) = getPhi(propPar.At(n, 3, 0), propPar.At(n, 4, 0));
    propPar_pol.At(n, 5, 0) = getTheta(pt, propPar.At(n, 5, 0));
    //jac_pol:first set to identity, then modify the relevant terms
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
	if (i==j) jac_pol.At(n, i, j) = 1.;
	else jac_pol.At(n, i, j) = 0.;
      }
    }
    jac_pol.At(n, 3, 3) = -propPar.At(n, 3, 0)/(pt*pt*pt);
    jac_pol.At(n, 3, 4) = -propPar.At(n, 4, 0)/(pt*pt*pt);
    jac_pol.At(n, 4, 3) = -propPar.At(n, 4, 0)/(pt*pt);
    jac_pol.At(n, 4, 4) =  propPar.At(n, 3, 0)/(pt*pt);
    float p2 = pt*pt + propPar.At(n, 5, 0)*propPar.At(n, 5, 0);
    jac_pol.At(n, 5, 3) =  propPar.At(n, 3, 0)*propPar.At(n, 5, 0)/(pt*p2);
    jac_pol.At(n, 5, 4) =  propPar.At(n, 4, 0)*propPar.At(n, 5, 0)/(pt*p2);
    jac_pol.At(n, 5, 5) = -pt/p2;
  }

#ifdef DEBUG
  if (dump) {
    printf("jac_pol:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", jac_pol.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexLS propErr_pol;
  MPlexLL propErr_tmp;
  MultFull(jac_pol,6,6,propErr,6,6,propErr_tmp,6,6);
  MultTranspFull(propErr_tmp,6,6,jac_pol,6,6,propErr_pol,6,6);

#ifdef DEBUG
  if (dump) {
    printf("propErr_pol:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", propErr_pol.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexLH K;
  MultFull(propErr_pol,6,6,resErrTmpLH,6,3,K,6,3);

#ifdef DEBUG
  if (dump) {
    printf("K:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", K.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexLV dPar;
  MultFull(K,6,3,res_loc,3,1,dPar,6,1);

  MPlexLV outPar_pol;
#pragma simd
  for (int n = 0; n < NN; ++n) {
    for (int i = 0; i < 6; ++i) {
      outPar_pol.At(n, i, 0) = propPar_pol.At(n, i, 0) + dPar.At(n, i, 0);
    }
  }

#ifdef DEBUG
  if (dump) {
    printf("outPar_pol:\n");
    for (int i = 0; i < 6; ++i) {
      printf("%8f  ", outPar_pol.At(0,i,0));
    } printf("\n");
    printf("\n");
  }
#endif

#pragma simd
  for (int n = 0; n < NN; ++n) {
    for (int i = 0; i < 3; ++i) {
      outPar.At(n, i, 0) = outPar_pol.At(n, i, 0);
    }
    outPar.At(n, 3, 0) = cos(outPar_pol.At(n, 4, 0))/outPar_pol.At(n, 3, 0);
    outPar.At(n, 4, 0) = sin(outPar_pol.At(n, 4, 0))/outPar_pol.At(n, 3, 0);
    outPar.At(n, 5, 0) = cos(outPar_pol.At(n, 5, 0))/(sin(outPar_pol.At(n, 5, 0))*outPar_pol.At(n, 3, 0));
  }

#ifdef DEBUG
  if (dump) {
    printf("outPar:\n");
    for (int i = 0; i < 6; ++i) {
      printf("%8f  ", outPar.At(0,i,0));
    } printf("\n");
    printf("\n");
  }
#endif

  MPlexLS I66;
#pragma simd
  for (int n = 0; n < NN; ++n) {
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
	if (i==j) I66.At(n, i, j) = 1.;
	else I66.At(n, i, j) = 0.;
      }
    }
  }

  MPlexHL H;
#pragma simd
  for (int n = 0; n < NN; ++n) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 6; ++j) {
	if (j>2) H.At(n, i, j) = 0.;
	else H.At(n, i, j) = rotT.At(n, i, j);
      }
    }
  }

#ifdef DEBUG
  if (dump) {
    printf("H:\n");
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", H.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexHS locErrMeas;
  MultFull      (rotT,3,3, msErr,3,3, temp,3,3);
  MultTranspFull(rotT,3,3, temp,3,3, locErrMeas,3,3);
#pragma simd
  for (int n = 0; n < NN; ++n) {
    locErrMeas.At(n, 2, 0) = 0;
    locErrMeas.At(n, 2, 1) = 0;
    locErrMeas.At(n, 2, 2) = 0;
    //locErrMeas.At(n, 1, 2) = 0;
    //locErrMeas.At(n, 0, 2) = 0;
  }

#ifdef DEBUG
  if (dump) {
    printf("locErrMeas:\n");
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j)
        printf("%8f ", locErrMeas.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexLS T1;
  MPlexLL tmp1;
  MPlexLL tmp11;
  MultFull(K,6,3,H,3,6,tmp1,6,6);
#pragma simd
  for (int n = 0; n < NN; ++n) {
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
	tmp1.At(n, i, j) = I66.At(n, i, j) - tmp1.At(n, i, j);
      }
    }
  }
  MultFull(tmp1,6,6,propErr_pol,6,6,tmp11,6,6);
  MultTranspFull(tmp11,6,6,tmp1,6,6,T1,6,6);

#ifdef DEBUG
  if (dump) {
    printf("T1:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", T1.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexLS T2;
  MPlexHL tmp2;
  MultTranspFull(locErrMeas,3,3,K,6,3,tmp2,3,6);

#ifdef DEBUG
  if (dump) {
    printf("tmp2:\n");
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", tmp2.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MultFull(K,6,3,tmp2,3,6,T2,6,6);

#ifdef DEBUG
  if (dump) {
    printf("T2:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", T2.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexLS outErr_pol;
#pragma simd
  for (int n = 0; n < NN; ++n) {
    for (int i = 0; i < 6; ++i) {
      for (int j = i; j < 6; ++j) {
	outErr_pol.At(n, i, j) = T1.At(n, i, j) + T2.At(n, i, j);
      }
    }
  }

#ifdef DEBUG
  if (dump) {
    printf("outErr_pol:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", outErr_pol.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexLL jac_back_pol;
#pragma simd
  for (int n = 0; n < NN; ++n) {
    //jac_back_pol:first set to identity, then modify the relevant terms
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
	if (i==j) jac_back_pol.At(n, i, j) = 1.;
	else jac_back_pol.At(n, i, j) = 0.;
      }
    }
    jac_back_pol.At(n, 3, 3) = -cos(outPar_pol.At(n, 4, 0))/pow(outPar_pol.At(n, 3, 0),2);
    jac_back_pol.At(n, 3, 4) = -sin(outPar_pol.At(n, 4, 0))/outPar_pol.At(n, 3, 0);
    jac_back_pol.At(n, 4, 3) = -sin(outPar_pol.At(n, 4, 0))/pow(outPar_pol.At(n, 3, 0),2);
    jac_back_pol.At(n, 4, 4) =  cos(outPar_pol.At(n, 4, 0))/outPar_pol.At(n, 3, 0);
    jac_back_pol.At(n, 5, 3) = -cos(outPar_pol.At(n, 5, 0))/(sin(outPar_pol.At(n, 5, 0))*pow(outPar_pol.At(n, 3, 0),2));
    jac_back_pol.At(n, 5, 5) = -1./(pow(sin(outPar_pol.At(n, 5, 0)),2)*outPar_pol.At(n, 3, 0));
  }

#ifdef DEBUG
  if (dump) {
    printf("jac_back_pol:\n");
    for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
        printf("%8f ", jac_back_pol.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  MPlexLL outErr_tmp;
  MultFull(jac_back_pol,6,6,outErr_pol,6,6,outErr_tmp,6,6);
  MultTranspFull(outErr_tmp,6,6,jac_back_pol,6,6,outErr,6,6);

#ifdef DEBUG
  if (dump) {
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

  // Also: resErr could be 3x3, kalmanGain 6x3

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

  MPlexHH rot;//we may not need this one
  MPlexHH rotT;
  MPlexHV res_glo;
  MPlexHS resErr_glo;
#pragma simd
  for (int n = 0; n < NN; ++n) {
    float r = hipo(msPar.ConstAt(n, 0, 0), msPar.ConstAt(n, 1, 0));
    rot.At(n, 0, 0) = -(msPar.ConstAt(n, 1, 0)+propPar.ConstAt(n, 1, 0))/(2*r);
    rot.At(n, 0, 1) = 0;
    rot.At(n, 0, 2) =  (msPar.ConstAt(n, 0, 0)+propPar.ConstAt(n, 0, 0))/(2*r);
    rot.At(n, 1, 0) = rot.ConstAt(n, 0, 2);
    rot.At(n, 1, 1) = 0;
    rot.At(n, 1, 2) = -rot.ConstAt(n, 0, 0);
    rot.At(n, 2, 0) = 0;
    rot.At(n, 2, 1) = 1;
    rot.At(n, 2, 2) = 0;
    //
    rotT.At(n, 0, 0) = rot.ConstAt(n, 0, 0);
    rotT.At(n, 0, 1) = rot.ConstAt(n, 1, 0);
    rotT.At(n, 0, 2) = rot.ConstAt(n, 2, 0);
    rotT.At(n, 1, 0) = rot.ConstAt(n, 0, 1);
    rotT.At(n, 1, 1) = rot.ConstAt(n, 1, 1);
    rotT.At(n, 1, 2) = rot.ConstAt(n, 2, 1);
    rotT.At(n, 2, 0) = rot.ConstAt(n, 0, 2);
    rotT.At(n, 2, 1) = rot.ConstAt(n, 1, 2);
    rotT.At(n, 2, 2) = rot.ConstAt(n, 2, 2);
    //
    for (int i = 0; i < 3; ++i) {
      res_glo.At(n, i, 0) = msPar.ConstAt(n, i, 0) - propPar.ConstAt(n, i, 0);
      for (int j = 0; j < 3; ++j) {
	resErr_glo.At(n, i, j) = msErr.ConstAt(n,i,j) + propErr.ConstAt(n,i,j);
      }
    }
  }

  MPlexHV res_loc;
  MultiplyGeneral(rotT,res_glo,res_loc);
  MPlexHS resErr_loc;
  MPlexHH temp;
  MultFull      (rotT,3,3, resErr_glo,3,3, temp,3,3);
  MultTranspFull(rotT,3,3, temp,3,3, resErr_loc,3,3);

  MPlex2S resErr;
#pragma simd
  for (int n = 0; n < NN; ++n) {
    for (int i = 0; i < 2; ++i) {
      for (int j = i; j < 2; ++j) {
	resErr.At(n, i, j) = resErr_loc.ConstAt(n,i,j);
      }
    }
  }

#ifdef DEBUG
  if (dump) {
    printf("resErr:\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  Matriplex::InvertCramerSym(resErr);

#ifdef DEBUG
  if (dump) {
    printf("resErrInv:\n");
    for (int i = 0; i < 2; ++i) { for (int j = 0; j < 2; ++j)
        printf("%8f ", resErr.At(0,i,j)); printf("\n");
    } printf("\n");
  }
#endif

  Chi2Similarity(res_loc, resErr, outChi2);
}
