#include "PropagationMPlex.h"

//==============================================================================
// propagateLineToRMPlex
//==============================================================================

using namespace Matriplex;

void propagateLineToRMPlex(const MPlexLS &psErr,  const MPlexLV& psPar,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                           MPlexLS &outErr,       MPlexLV& outPar)
{
   // XXX Regenerate parts below with a script.

   const idx_t N  = NN;

#pragma simd
   for (int n = 0; n < N; ++n)
   {

     const float cosA = (psPar[0 * N + n] * psPar[3 * N + n] + psPar[1 * N + n] * psPar[4 * N + n]) / ( sqrt( ( psPar[0 * N + n] * psPar[0 * N + n] + psPar[1 * N + n] * psPar[1 * N + n] ) * ( psPar[3 * N + n] * psPar[3 * N + n] + psPar[4 * N + n] * psPar[4 * N + n] ) ) );
     const float dr  = (hipo(msPar[0 * N + n], msPar[1 * N + n]) - hipo(psPar[0 * N + n], psPar[1 * N + n])) / cosA;

#ifdef DEBUG
     std::cout << "propagateLineToRMPlex dr=" << dr << std::endl;
#endif

      const float pt  = hipo(psPar[3 * N + n], psPar[4 * N + n]);
      const float p   = dr / pt; // path
      const float psq = p * p;

      outPar[0 * N + n] = psPar[0 * N + n] + p * psPar[3 * N + n];
      outPar[1 * N + n] = psPar[1 * N + n] + p * psPar[4 * N + n];
      outPar[2 * N + n] = psPar[2 * N + n] + p * psPar[5 * N + n];
      outPar[3 * N + n] = psPar[3 * N + n];
      outPar[4 * N + n] = psPar[4 * N + n];
      outPar[5 * N + n] = psPar[5 * N + n];

      {
        const MPlexLS& A = psErr;
              MPlexLS& B = outErr;

        B.fArray[0 * N + n] = A.fArray[0 * N + n];
        B.fArray[1 * N + n] = A.fArray[1 * N + n];
        B.fArray[2 * N + n] = A.fArray[2 * N + n];
        B.fArray[3 * N + n] = A.fArray[3 * N + n];
        B.fArray[4 * N + n] = A.fArray[4 * N + n];
        B.fArray[5 * N + n] = A.fArray[5 * N + n];
        B.fArray[6 * N + n] = A.fArray[6 * N + n] + p * A.fArray[0 * N + n];
        B.fArray[7 * N + n] = A.fArray[7 * N + n] + p * A.fArray[1 * N + n];
        B.fArray[8 * N + n] = A.fArray[8 * N + n] + p * A.fArray[3 * N + n];
        B.fArray[9 * N + n] = A.fArray[9 * N + n] + p * (A.fArray[6 * N + n] + A.fArray[6 * N + n]) + psq * A.fArray[0 * N + n];
        B.fArray[10 * N + n] = A.fArray[10 * N + n] + p * A.fArray[1 * N + n];
        B.fArray[11 * N + n] = A.fArray[11 * N + n] + p * A.fArray[2 * N + n];
        B.fArray[12 * N + n] = A.fArray[12 * N + n] + p * A.fArray[4 * N + n];
        B.fArray[13 * N + n] = A.fArray[13 * N + n] + p * (A.fArray[7 * N + n] + A.fArray[10 * N + n]) + psq * A.fArray[1 * N + n];
        B.fArray[14 * N + n] = A.fArray[14 * N + n] + p * (A.fArray[11 * N + n] + A.fArray[11 * N + n]) + psq * A.fArray[2 * N + n];
        B.fArray[15 * N + n] = A.fArray[15 * N + n] + p * A.fArray[3 * N + n];
        B.fArray[16 * N + n] = A.fArray[16 * N + n] + p * A.fArray[4 * N + n];
        B.fArray[17 * N + n] = A.fArray[17 * N + n] + p * A.fArray[5 * N + n];
        B.fArray[18 * N + n] = A.fArray[18 * N + n] + p * (A.fArray[8 * N + n] + A.fArray[15 * N + n]) + psq * A.fArray[3 * N + n];
        B.fArray[19 * N + n] = A.fArray[19 * N + n] + p * (A.fArray[12 * N + n] + A.fArray[16 * N + n]) + psq * A.fArray[4 * N + n];
        B.fArray[20 * N + n] = A.fArray[20 * N + n] + p * (A.fArray[17 * N + n] + A.fArray[17 * N + n]) + psq * A.fArray[5 * N + n];
      }

#ifdef DEBUG
      std::cout << "propagateLineToRMPlex arrive at r=" << hipo(outPar[0 * N + n], outPar[1 * N + n]) << std::endl;
#endif

   }
}


//==============================================================================
// propagateHelixToRMPlex
//==============================================================================

namespace
{

void MultHelixProp(const MPlexLL& A, const MPlexLS& B, MPlexLL& C)
{
   // C = A * B

   typedef float T;
   const idx_t N  = NN;

   const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
   const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
         T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "MultHelixProp.ah"
}

void MultHelixPropTransp(const MPlexLL& A, const MPlexLL& B, MPlexLS& C)
{
   // C = B * AT;

   typedef float T;
   const idx_t N  = NN;

   const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
   const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
         T *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "MultHelixPropTransp.ah"
}

inline void MultHelixPropTemp(const MPlexLL& A, const MPlexLL& B, MPlexLL& C, int n)
{
   // C = A * B

  typedef float T;
  const idx_t N  = NN;

  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const T *b = B.fArray; ASSUME_ALIGNED(b, 64);
  T *c = C.fArray; ASSUME_ALIGNED(c, 64);

  c[ 0*N+n] = a[ 0*N+n]*b[ 0*N+n] + a[ 1*N+n]*b[ 6*N+n] + a[ 2*N+n]*b[12*N+n] + a[ 4*N+n]*b[24*N+n];
  c[ 1*N+n] = a[ 0*N+n]*b[ 1*N+n] + a[ 1*N+n]*b[ 7*N+n] + a[ 2*N+n]*b[13*N+n] + a[ 4*N+n]*b[25*N+n];
  c[ 2*N+n] = a[ 2*N+n];
  c[ 3*N+n] = a[ 0*N+n]*b[ 3*N+n] + a[ 1*N+n]*b[ 9*N+n] + a[ 2*N+n]*b[15*N+n] + a[ 3*N+n] + a[ 4*N+n]*b[27*N+n];
  c[ 4*N+n] = a[ 0*N+n]*b[ 4*N+n] + a[ 1*N+n]*b[10*N+n] + a[ 4*N+n];
  c[ 5*N+n] = a[ 2*N+n]*b[17*N+n] + a[ 5*N+n];
  c[ 6*N+n] = a[ 6*N+n]*b[ 0*N+n] + a[ 7*N+n]*b[ 6*N+n] + a[ 8*N+n]*b[12*N+n] + a[10*N+n]*b[24*N+n];
  c[ 7*N+n] = a[ 6*N+n]*b[ 1*N+n] + a[ 7*N+n]*b[ 7*N+n] + a[ 8*N+n]*b[13*N+n] + a[10*N+n]*b[25*N+n];
  c[ 8*N+n] = a[ 8*N+n];
  c[ 9*N+n] = a[ 6*N+n]*b[ 3*N+n] + a[ 7*N+n]*b[ 9*N+n] + a[ 8*N+n]*b[15*N+n] + a[ 9*N+n] + a[10*N+n]*b[27*N+n];
  c[10*N+n] = a[ 6*N+n]*b[ 4*N+n] + a[ 7*N+n]*b[10*N+n] + a[10*N+n];
  c[11*N+n] = a[ 8*N+n]*b[17*N+n] + a[11*N+n];
  c[12*N+n] = a[12*N+n]*b[ 0*N+n] + a[13*N+n]*b[ 6*N+n] + a[14*N+n]*b[12*N+n] + a[16*N+n]*b[24*N+n];
  c[13*N+n] = a[12*N+n]*b[ 1*N+n] + a[13*N+n]*b[ 7*N+n] + a[14*N+n]*b[13*N+n] + a[16*N+n]*b[25*N+n];
  c[14*N+n] = a[14*N+n];
  c[15*N+n] = a[12*N+n]*b[ 3*N+n] + a[13*N+n]*b[ 9*N+n] + a[14*N+n]*b[15*N+n] + a[15*N+n] + a[16*N+n]*b[27*N+n];
  c[16*N+n] = a[12*N+n]*b[ 4*N+n] + a[13*N+n]*b[10*N+n] + a[16*N+n];
  c[17*N+n] = a[14*N+n]*b[17*N+n] + a[17*N+n];
  c[18*N+n] = a[18*N+n]*b[ 0*N+n] + a[19*N+n]*b[ 6*N+n] + a[20*N+n]*b[12*N+n] + a[22*N+n]*b[24*N+n];
  c[19*N+n] = a[18*N+n]*b[ 1*N+n] + a[19*N+n]*b[ 7*N+n] + a[20*N+n]*b[13*N+n] + a[22*N+n]*b[25*N+n];
  c[20*N+n] = a[20*N+n];
  c[21*N+n] = a[18*N+n]*b[ 3*N+n] + a[19*N+n]*b[ 9*N+n] + a[20*N+n]*b[15*N+n] + a[21*N+n] + a[22*N+n]*b[27*N+n];
  c[22*N+n] = a[18*N+n]*b[ 4*N+n] + a[19*N+n]*b[10*N+n] + a[22*N+n];
  c[23*N+n] = a[20*N+n]*b[17*N+n] + a[23*N+n];
  c[24*N+n] = a[24*N+n]*b[ 0*N+n] + a[25*N+n]*b[ 6*N+n] + a[26*N+n]*b[12*N+n] + a[28*N+n]*b[24*N+n];
  c[25*N+n] = a[24*N+n]*b[ 1*N+n] + a[25*N+n]*b[ 7*N+n] + a[26*N+n]*b[13*N+n] + a[28*N+n]*b[25*N+n];
  c[26*N+n] = a[26*N+n];
  c[27*N+n] = a[24*N+n]*b[ 3*N+n] + a[25*N+n]*b[ 9*N+n] + a[26*N+n]*b[15*N+n] + a[27*N+n] + a[28*N+n]*b[27*N+n];
  c[28*N+n] = a[24*N+n]*b[ 4*N+n] + a[25*N+n]*b[10*N+n] + a[28*N+n];
  c[29*N+n] = a[26*N+n]*b[17*N+n] + a[29*N+n];
  c[30*N+n] = a[30*N+n]*b[ 0*N+n] + a[31*N+n]*b[ 6*N+n] + a[32*N+n]*b[12*N+n] + a[34*N+n]*b[24*N+n];
  c[31*N+n] = a[30*N+n]*b[ 1*N+n] + a[31*N+n]*b[ 7*N+n] + a[32*N+n]*b[13*N+n] + a[34*N+n]*b[25*N+n];
  c[32*N+n] = a[32*N+n];
  c[33*N+n] = a[30*N+n]*b[ 3*N+n] + a[31*N+n]*b[ 9*N+n] + a[32*N+n]*b[15*N+n] + a[33*N+n] + a[34*N+n]*b[27*N+n];
  c[34*N+n] = a[30*N+n]*b[ 4*N+n] + a[31*N+n]*b[10*N+n] + a[34*N+n];
  c[35*N+n] = a[32*N+n]*b[17*N+n] + a[35*N+n];
}

// this version does not assume to know which elements are 0 or 1, so it does the full multiplication
void MultHelixPropFull(const MPlexLL& A, const MPlexLS& B, MPlexLL& C)
{

#pragma simd
  for (int n = 0; n < NN; ++n)
    {
      for (int i = 0; i < 6; ++i) {
	for (int j = 0; j < 6; ++j) {
	  C(n,i,j) = 0.;
	  for (int k = 0; k < 6; ++k) C(n,i,j) += A.ConstAt(n,i,k)*B.ConstAt(n,k,j);
	}
      }
    }
  
}

// this version does not assume to know which elements are 0 or 1, so it does the full multiplication
void MultHelixPropFull(const MPlexLL& A, const MPlexLL& B, MPlexLL& C)
{

#pragma simd
  for (int n = 0; n < NN; ++n)
    {
      for (int i = 0; i < 6; ++i) {
	for (int j = 0; j < 6; ++j) {
	  C(n,i,j) = 0.;
	  for (int k = 0; k < 6; ++k) C(n,i,j) += A.ConstAt(n,i,k)*B.ConstAt(n,k,j);
	}
      }
    }
  
}

// this version does not assume to know which elements are 0 or 1, so it does the full mupltiplication
void MultHelixPropTranspFull(const MPlexLL& A, const MPlexLL& B, MPlexLS& C)
{

#pragma simd
  for (int n = 0; n < NN; ++n)
    {
      for (int i = 0; i < 6; ++i) {
	for (int j = 0; j < 6; ++j) {
	  C(n,i,j) = 0.;
	  for (int k = 0; k < 6; ++k) C(n,i,j) += B.ConstAt(n,i,k)*A.ConstAt(n,j,k);
	}
      }
    }
  
}

// this version does not assume to know which elements are 0 or 1, so it does the full mupltiplication
void MultHelixPropTranspFull(const MPlexLL& A, const MPlexLL& B, MPlexLL& C)
{

#pragma simd
  for (int n = 0; n < NN; ++n)
    {
      for (int i = 0; i < 6; ++i) {
	for (int j = 0; j < 6; ++j) {
	  C(n,i,j) = 0.;
	  for (int k = 0; k < 6; ++k) C(n,i,j) += B.ConstAt(n,i,k)*A.ConstAt(n,j,k);
	}
      }
    }
  
}

}

//#define DEBUG

void helixAtRFromIterativePolarFullJac(const MPlexLV& inPar, const MPlexQI& inChg, MPlexLV& outPar, const MPlexQF &msRad, MPlexLL& errorProp) {

  errorProp.SetVal(0);
  MPlexLL errorPropTmp(0);//initialize to zero
  MPlexLL errorPropSwap(0);//initialize to zero

#pragma simd
  for (int n = 0; n < NN; ++n)
    {

      //initialize erroProp to identity matrix
      errorProp(n,0,0) = 1.;
      errorProp(n,1,1) = 1.;
      errorProp(n,2,2) = 1.;
      errorProp(n,3,3) = 1.;
      errorProp(n,4,4) = 1.;
      errorProp(n,5,5) = 1.;

      const float k  = inChg.ConstAt(n, 0, 0) * 100. / (-Config::sol*Config::Bfield);
      const float& r = msRad.ConstAt(n, 0, 0);
      float r0 = hipo(inPar.ConstAt(n, 0, 0), inPar.ConstAt(n, 1, 0));

      if (fabs(r-r0)<0.0001) {
#ifdef DEBUG
	std::cout << "distance less than 1mum, skip" << std::endl;
#endif
	continue;
      }

      const float& xin   = inPar.ConstAt(n, 0, 0);
      const float& yin   = inPar.ConstAt(n, 1, 0);
      const float& zin   = inPar.ConstAt(n, 2, 0);
      const float& ipt   = inPar.ConstAt(n, 3, 0);
      const float& phiin = inPar.ConstAt(n, 4, 0);
      const float& theta = inPar.ConstAt(n, 5, 0);

      //set those that are 1. before iterations
      errorPropTmp(n,2,2) = 1.;
      errorPropTmp(n,3,3) = 1.;
      errorPropTmp(n,4,4) = 1.;
      errorPropTmp(n,5,5) = 1.;

      float /*alpha  = 0.,*/ cosa = 0., sina = 0., ialpha = 0.;
      float cosP = 0., sinP = 0.;
      float cosT = 0., sinT = 0.;
      if (Config::useTrigApprox) {
	sincos4(phiin, sinP, cosP);
	sincos4(theta, sinT, cosT);
      } else {
	cosP=cos(phiin);
	sinP=sin(phiin);
	cosT=cos(theta);
	sinT=sin(theta);
      }
      float pxin = cosP/ipt;
      float pyin = sinP/ipt;
      float dadx = 0., dady = 0., dadipt = 0.;

      for (unsigned int i=0;i<Config::Niter;++i) {

#ifdef DEBUG
	std::cout << std::endl;
	std::cout << "attempt propagation from r=" << r0 << " to r=" << r << std::endl;
	std::cout << "x=" << outPar.At(n, 0, 0) << " y=" << outPar.At(n, 1, 0)  << " z=" << outPar.At(n, 2, 0)
		  << " px=" << cos(phiin)/ipt << " py=" << sin(phiin)/ipt << " pz=" << 1./(ipt*tan(theta)) << " q=" << inChg.ConstAt(n, 0, 0) << std::endl;
#endif

	r0 = hipo(outPar.ConstAt(n, 0, 0), outPar.ConstAt(n, 1, 0));
	ialpha = (r-r0)*ipt/k;
	//alpha+=ialpha;

	if (Config::useTrigApprox) {
	  sincos4(ialpha, sina, cosa);
	} else {
	  cosa=cos(ialpha);
	  sina=sin(ialpha);
	}

	//derivatives of alpha
	dadx   = -outPar.At(n, 0, 0)*ipt/(k*r0);
	dady   = -outPar.At(n, 1, 0)*ipt/(k*r0);
	dadipt = (r-r0)/k;

	outPar.At(n, 0, 0) = outPar.ConstAt(n, 0, 0) + k*(pxin*sina-pyin*(1.-cosa));
	outPar.At(n, 1, 0) = outPar.ConstAt(n, 1, 0) + k*(pyin*sina+pxin*(1.-cosa));
	const float pxinold = pxin;//copy before overwriting
	pxin = pxin*cosa-pyin*sina;
	pyin = pyin*cosa+pxinold*sina;

	//need phi at origin, so this goes before redefining phi
	if (Config::useTrigApprox) {
	  sincos4(outPar.At(n, 4, 0), sinP, cosP);
	} else {
	  cosP=cos(outPar.At(n, 4, 0));
	  sinP=sin(outPar.At(n, 4, 0));
	}

	outPar.At(n, 2, 0) = outPar.ConstAt(n, 2, 0) + k*ialpha*cosT/(ipt*sinT);
	outPar.At(n, 3, 0) = ipt;
	outPar.At(n, 4, 0) = outPar.ConstAt(n, 4, 0)+ialpha;
	outPar.At(n, 5, 0) = theta;

	errorPropTmp(n,0,0) = 1.+k*(cosP*dadx*cosa-sinP*dadx*sina)/ipt;
	errorPropTmp(n,0,1) =    k*(cosP*dady*cosa-sinP*dady*sina)/ipt;
	errorPropTmp(n,0,3) = k*(cosP*(ipt*dadipt*cosa-sina)+sinP*((1.-cosa)-ipt*dadipt*sina))/(ipt*ipt);
	errorPropTmp(n,0,4) =-k*(sinP*sina+cosP*(1.-cosa))/ipt;

	errorPropTmp(n,1,0) =    k*(sinP*dadx*cosa+cosP*dadx*sina)/ipt;
	errorPropTmp(n,1,1) = 1.+k*(sinP*dady*cosa+cosP*dady*sina)/ipt;
	errorPropTmp(n,1,3) = k*(sinP*(ipt*dadipt*cosa-sina)+cosP*(ipt*dadipt*sina-(1.-cosa)))/(ipt*ipt);
	errorPropTmp(n,1,4) = k*(cosP*sina-sinP*(1.-cosa))/ipt;

	errorPropTmp(n,2,0) = k*cosT*dadx/(ipt*sinT);
	errorPropTmp(n,2,1) = k*cosT*dady/(ipt*sinT);
	errorPropTmp(n,2,3) = k*cosT*(ipt*dadipt-ialpha)/(ipt*ipt*sinT);
	errorPropTmp(n,2,5) =-k*ialpha/(ipt*sinT*sinT);

	errorPropTmp(n,4,0) = dadx;
	errorPropTmp(n,4,1) = dady;
	errorPropTmp(n,4,3) = dadipt;

	MultHelixPropTemp(errorProp,errorPropTmp,errorPropSwap,n);
	errorProp = errorPropSwap;

      }

#ifdef DEBUG
      std::cout << "propagation end, dump parameters" << std::endl;
      std::cout << "pos = " << outPar.At(n, 0, 0) << " " << outPar.At(n, 1, 0) << " " << outPar.At(n, 2, 0) << std::endl;
      std::cout << "mom = " << cos(outPar.At(n, 4, 0))/outPar.At(n, 3, 0) << " " << sin(outPar.At(n, 4, 0))/outPar.At(n, 3, 0) << " " << 1./(outPar.At(n, 3, 0)*tan(outPar.At(n, 5, 0))) << std::endl;
      std::cout << "r=" << sqrt( outPar.At(n, 0, 0)*outPar.At(n, 0, 0) + outPar.At(n, 1, 0)*outPar.At(n, 1, 0) ) << " pT=" << 1./fabs(outPar.At(n, 3, 0)) << std::endl;

      std::cout << "jacobian" << std::endl;
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,0,0),errorProp(n,0,1),errorProp(n,0,2),errorProp(n,0,3),errorProp(n,0,4),errorProp(n,0,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,1,0),errorProp(n,1,1),errorProp(n,1,2),errorProp(n,1,3),errorProp(n,1,4),errorProp(n,1,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,2,0),errorProp(n,2,1),errorProp(n,2,2),errorProp(n,2,3),errorProp(n,2,4),errorProp(n,2,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,3,0),errorProp(n,3,1),errorProp(n,3,2),errorProp(n,3,3),errorProp(n,3,4),errorProp(n,3,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,4,0),errorProp(n,4,1),errorProp(n,4,2),errorProp(n,4,3),errorProp(n,4,4),errorProp(n,4,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,5,0),errorProp(n,5,1),errorProp(n,5,2),errorProp(n,5,3),errorProp(n,5,4),errorProp(n,5,5));
#endif

    }
}

void helixAtRFromIterativePolar(const MPlexLV& inPar, const MPlexQI& inChg, MPlexLV& outPar, const MPlexQF &msRad, MPlexLL& errorProp) {

  errorProp.SetVal(0);

#pragma simd
  for (int n = 0; n < NN; ++n)
    {

      //initialize erroProp to identity matrix
      errorProp(n,0,0) = 1.;
      errorProp(n,1,1) = 1.;
      errorProp(n,2,2) = 1.;
      errorProp(n,3,3) = 1.;
      errorProp(n,4,4) = 1.;
      errorProp(n,5,5) = 1.;

      const float k  = inChg.ConstAt(n, 0, 0) * 100. / (-Config::sol*Config::Bfield);
      const float& r = msRad.ConstAt(n, 0, 0);
      float r0 = hipo(inPar.ConstAt(n, 0, 0), inPar.ConstAt(n, 1, 0));

      if (fabs(r-r0)<0.0001) {
#ifdef DEBUG
	std::cout << "distance less than 1mum, skip" << std::endl;
#endif
	continue;
      }

      const float& xin   = inPar.ConstAt(n, 0, 0);
      const float& yin   = inPar.ConstAt(n, 1, 0);
      const float& zin   = inPar.ConstAt(n, 2, 0);
      const float& ipt   = inPar.ConstAt(n, 3, 0);
      const float& phiin = inPar.ConstAt(n, 4, 0);
      const float& theta = inPar.ConstAt(n, 5, 0);

      const float kinv  = 1.f/k;
      const float pt = 1.f/ipt;

      float D = 0., cosa = 0., sina = 0., id = 0.;
      float cosPorT = 0., sinPorT = 0.;
      if (Config::useTrigApprox) {
	sincos4(phiin, sinPorT, cosPorT);
      } else {
	cosPorT=cos(phiin);
	sinPorT=sin(phiin);
      }
      float pxin = cosPorT*pt;
      float pyin = sinPorT*pt;

      //derivatives initialized to value for first iteration, i.e. distance = r-r0in
      float dDdx = r0>0. ? -xin/r0 : 0.;
      float dDdy = r0>0. ? -yin/r0 : 0.;
      float dDdipt = 0.;
      float dDdphi = 0.;
      //temporaries
      float dadx = 0.;
      float dady = 0.;
      float dadipt = 0.;

      for (unsigned int i=0;i<Config::Niter;++i) {

#ifdef DEBUG
	std::cout << std::endl;
	std::cout << "attempt propagation from r=" << r0 << " to r=" << r << std::endl;
	std::cout << "x=" << outPar.At(n, 0, 0) << " y=" << outPar.At(n, 1, 0)  << " z=" << outPar.At(n, 2, 0)
		  << " px=" << cos(phiin)*pt << " py=" << sin(phiin)*pt << " pz=" << 1./(ipt*tan(theta)) << " q=" << inChg.ConstAt(n, 0, 0) << std::endl;
#endif

	//compute distance and path for the current iteration
	r0 = hipo(outPar.At(n, 0, 0), outPar.At(n, 1, 0));
	id = (r-r0);
	D+=id;
	if (Config::useTrigApprox) {
	  sincos4(id*ipt*kinv, sina, cosa);
	} else {
	  cosa=cos(id*ipt*kinv);
	  sina=sin(id*ipt*kinv);
	}

	//update derivatives on total distance
	if (i+1 != Config::Niter && r0 > 0 && fabs((r-r0))>0.0001) {

	  float& x = outPar.At(n, 0, 0);
	  float& y = outPar.At(n, 1, 0);

	  //redefine r0 as 1./r0 to reduce the number of temporaries
	  r0 = 1./r0;

	  dadx = -x*ipt*kinv*r0;
	  dady = -y*ipt*kinv*r0;
	  dadipt = id*kinv;

	  dDdx   -= ( x*(1.f + k*dadx*(pxin*cosa - pyin*sina)) + y*k*dadx*(pyin*cosa + pxin*sina) )*r0;
	  dDdy   -= ( x*k*dady*(pxin*cosa - pyin*sina) + y*(1.f + k*dady*(pyin*cosa + pxin*sina)) )*r0;
	  //now r0 depends on ipt and phi as well
	  dDdipt -= k*( x*(pxin*dadipt*cosa*ipt - pyin*dadipt*sina*ipt - pyin*cosa - pxin*sina + pyin) +
			y*(pyin*dadipt*cosa*ipt + pxin*dadipt*sina*ipt - pyin*sina + pxin*cosa - pxin))*pt*r0;
	  dDdphi += k*( x*(pyin*sina - pxin + pxin*cosa) - y*(pxin*sina - pyin + pyin*cosa))*r0;

	}

	//update parameters
	outPar.At(n, 0, 0) = outPar.At(n, 0, 0) + k*(pxin*sina - pyin*(1.-cosa));
	outPar.At(n, 1, 0) = outPar.At(n, 1, 0) + k*(pyin*sina + pxin*(1.-cosa));
	const float pxinold = pxin;//copy before overwriting
	pxin = pxin*cosa - pyin*sina;
	pyin = pyin*cosa + pxinold*sina;

      }

      const float alpha  = D*ipt*kinv;
      dadx   = dDdx*ipt*kinv;
      dady   = dDdy*ipt*kinv;
      dadipt = (ipt*dDdipt + D)*kinv;
      const float dadphi = dDdphi*ipt*kinv;

      if (Config::useTrigApprox) {
	sincos4(alpha, sina, cosa);
      } else {
	cosa=cos(alpha);
	sina=sin(alpha);
      }

      errorProp(n,0,0) = 1.+k*dadx*(cosPorT*cosa-sinPorT*sina)*pt;
      errorProp(n,0,1) =    k*dady*(cosPorT*cosa-sinPorT*sina)*pt;
      errorProp(n,0,2) = 0.;
      errorProp(n,0,3) = k*(cosPorT*(ipt*dadipt*cosa-sina)+sinPorT*((1.-cosa)-ipt*dadipt*sina))*pt*pt;
      errorProp(n,0,4) = k*(cosPorT*dadphi*cosa - sinPorT*dadphi*sina - sinPorT*sina + cosPorT*cosa - cosPorT)*pt;
      errorProp(n,0,5) = 0.;

      errorProp(n,1,0) =    k*dadx*(sinPorT*cosa+cosPorT*sina)*pt;
      errorProp(n,1,1) = 1.+k*dady*(sinPorT*cosa+cosPorT*sina)*pt;
      errorProp(n,1,2) = 0.;
      errorProp(n,1,3) = k*(sinPorT*(ipt*dadipt*cosa-sina)+cosPorT*(ipt*dadipt*sina-(1.-cosa)))*pt*pt;
      errorProp(n,1,4) = k*(sinPorT*dadphi*cosa + cosPorT*dadphi*sina + sinPorT*cosa + cosPorT*sina - sinPorT)*pt;
      errorProp(n,1,5) = 0.;

      if (Config::useTrigApprox) {
	sincos4(theta, sinPorT, cosPorT);
      } else {
	cosPorT=cos(theta);
	sinPorT=sin(theta);
      }
      //redefine sinPorT as 1./sinPorT to reduce the number of temporaries
      sinPorT = 1./sinPorT;

      outPar.At(n, 2, 0) = inPar.ConstAt(n, 2, 0) + k*alpha*cosPorT*pt*sinPorT;

      errorProp(n,2,0) = k*cosPorT*dadx*pt*sinPorT;
      errorProp(n,2,1) = k*cosPorT*dady*pt*sinPorT;
      errorProp(n,2,2) = 1.;
      errorProp(n,2,3) = k*cosPorT*(ipt*dadipt-alpha)*pt*pt*sinPorT;
      errorProp(n,2,4) = k*dadphi*cosPorT*pt*sinPorT;
      errorProp(n,2,5) =-k*alpha*pt*sinPorT*sinPorT;

      outPar.At(n, 3, 0) = ipt;

      errorProp(n,3,0) = 0.;
      errorProp(n,3,1) = 0.;
      errorProp(n,3,2) = 0.;
      errorProp(n,3,3) = 1.;
      errorProp(n,3,4) = 0.;
      errorProp(n,3,5) = 0.;

      outPar.At(n, 4, 0) = inPar.ConstAt(n, 4, 0)+alpha;

      errorProp(n,4,0) = dadx;
      errorProp(n,4,1) = dady;
      errorProp(n,4,2) = 0.;
      errorProp(n,4,3) = dadipt;
      errorProp(n,4,4) = 1.+dadphi;
      errorProp(n,4,5) = 0.;

      outPar.At(n, 5, 0) = theta;

      errorProp(n,5,0) = 0.;
      errorProp(n,5,1) = 0.;
      errorProp(n,5,2) = 0.;
      errorProp(n,5,3) = 0.;
      errorProp(n,5,4) = 0.;
      errorProp(n,5,5) = 1.;

#ifdef DEBUG
      std::cout << "propagation end, dump parameters" << std::endl;
      std::cout << "pos = " << outPar.At(n, 0, 0) << " " << outPar.At(n, 1, 0) << " " << outPar.At(n, 2, 0) << std::endl;
      std::cout << "mom = " << cos(outPar.At(n, 4, 0))/outPar.At(n, 3, 0) << " " << sin(outPar.At(n, 4, 0))/outPar.At(n, 3, 0) << " " << 1./(outPar.At(n, 3, 0)*tan(outPar.At(n, 5, 0))) << std::endl;
      std::cout << "r=" << sqrt( outPar.At(n, 0, 0)*outPar.At(n, 0, 0) + outPar.At(n, 1, 0)*outPar.At(n, 1, 0) ) << " pT=" << 1./fabs(outPar.At(n, 3, 0)) << std::endl;

      std::cout << "jacobian" << std::endl;
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,0,0),errorProp(n,0,1),errorProp(n,0,2),errorProp(n,0,3),errorProp(n,0,4),errorProp(n,0,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,1,0),errorProp(n,1,1),errorProp(n,1,2),errorProp(n,1,3),errorProp(n,1,4),errorProp(n,1,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,2,0),errorProp(n,2,1),errorProp(n,2,2),errorProp(n,2,3),errorProp(n,2,4),errorProp(n,2,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,3,0),errorProp(n,3,1),errorProp(n,3,2),errorProp(n,3,3),errorProp(n,3,4),errorProp(n,3,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,4,0),errorProp(n,4,1),errorProp(n,4,2),errorProp(n,4,3),errorProp(n,4,4),errorProp(n,4,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,5,0),errorProp(n,5,1),errorProp(n,5,2),errorProp(n,5,3),errorProp(n,5,4),errorProp(n,5,5));
#endif

    }
}

void helixAtRFromIterative(const MPlexLV& inPar, const MPlexQI& inChg, MPlexLV& outPar, const MPlexQF &msRad, MPlexLL& errorProp) {

  errorProp.SetVal(0);

#pragma simd
  for (int n = 0; n < NN; ++n)
    {

      //initialize erroProp to identity matrix
      errorProp(n,0,0) = 1.;
      errorProp(n,1,1) = 1.;
      errorProp(n,2,2) = 1.;
      errorProp(n,3,3) = 1.;
      errorProp(n,4,4) = 1.;
      errorProp(n,5,5) = 1.;

      const float& xin  = inPar.ConstAt(n, 0, 0);
      const float& yin  = inPar.ConstAt(n, 1, 0);
      const float& pxin = inPar.ConstAt(n, 3, 0);
      const float& pyin = inPar.ConstAt(n, 4, 0);
      const float& pzin = inPar.ConstAt(n, 5, 0);
      const float& r    = msRad.ConstAt(n, 0, 0);
      float r0 = hipo(xin, yin);
      
#ifdef DEBUG
      std::cout << std::endl;
      std::cout << "attempt propagation from r=" << r0 << " to r=" << r << std::endl;
      std::cout << "x=" << xin << " y=" << yin  << " z=" << inPar.ConstAt(n, 2, 0) << " px=" << pxin << " py=" << pyin << " pz=" << pzin << " q=" << inChg.ConstAt(n, 0, 0) << std::endl;
#endif

      if (fabs(r-r0)<0.0001) {
#ifdef DEBUG
	std::cout << "distance less than 1mum, skip" << std::endl;
#endif
	continue;
      }
      
      float pt2    = pxin*pxin+pyin*pyin;
      float pt     = sqrt(pt2);
      float ptinv  = 1./pt;
      float pt2inv = ptinv*ptinv;
      //p=0.3Br => r=p/(0.3*B)
      float k = inChg.ConstAt(n, 0, 0) * 100. / (-Config::sol*Config::Bfield);
      float invcurvature = 1./(pt*k);//in 1./cm
      float ctgTheta=pzin*ptinv;
      
#ifdef DEBUG
      std::cout << "curvature=" << 1./invcurvature << std::endl;
#endif
      
      //variables to be updated at each iterations
      float totalDistance = 0;
      //derivatives initialized to value for first iteration, i.e. distance = r-r0in
      float dTDdx = r0>0. ? -xin/r0 : 0.;
      float dTDdy = r0>0. ? -yin/r0 : 0.;
      float dTDdpx = 0.;
      float dTDdpy = 0.;
      //temporaries used within the loop (declare here to reduce memory operations)
      float x = 0.;
      float y = 0.;
      float px = 0.;
      float py = 0.;
      float cosAP=0.;
      float sinAP=0.;
      float dAPdx = 0.;
      float dAPdy = 0.;
      float dAPdpx = 0.;
      float dAPdpy = 0.;
      for (unsigned int i=0;i<Config::Niter;++i)
	{
#ifdef DEBUG
	  std::cout << "propagation iteration #" << i << std::endl;
#endif
	  
	  x  = outPar.At(n, 0, 0);
	  y  = outPar.At(n, 1, 0);
	  px = outPar.At(n, 3, 0);
	  py = outPar.At(n, 4, 0);
	  r0 = hipo(outPar.At(n, 0, 0), outPar.At(n, 1, 0));
	  
#ifdef DEBUG
	  std::cout << "r0=" << r0 << " pt=" << pt << std::endl;
	  // if (dump) {
	  //    if (r==r0) {
	  //       std::cout << "distance = 0 at iteration=" << i << std::endl;
	  //       break;
	  //    }
	  // }
#endif
	  
	  //distance=r-r0;//remove temporary
	  totalDistance+=(r-r0);
	  
#ifdef DEBUG
	  std::cout << "distance=" << (r-r0) << " angPath=" << (r-r0)*invcurvature << std::endl;
#endif
	  
	  //float angPath = (r-r0)*invcurvature;
	  if (Config::useTrigApprox) {
	    sincos4((r-r0)*invcurvature, sinAP, cosAP);
	  } else {
	    cosAP=cos((r-r0)*invcurvature);
	    sinAP=sin((r-r0)*invcurvature);
	  }

	  //helix propagation formulas
	  //http://www.phys.ufl.edu/~avery/fitting/fitting4.pdf
	  outPar.At(n, 0, 0) = outPar.At(n, 0, 0) + k*(px*sinAP-py*(1-cosAP));
	  outPar.At(n, 1, 0) = outPar.At(n, 1, 0) + k*(py*sinAP+px*(1-cosAP));
	  outPar.At(n, 2, 0) = outPar.At(n, 2, 0) + (r-r0)*ctgTheta;
	  outPar.At(n, 3, 0) = px*cosAP-py*sinAP;
	  outPar.At(n, 4, 0) = py*cosAP+px*sinAP;
	  //outPar.At(n, 5, 0) = pz; //take this out as it is redundant

	  if (i+1 != Config::Niter && r0 > 0 && fabs((r-r0)*invcurvature)>0.000000001)
	  {
	     //update derivatives on total distance for next step, where totalDistance+=r-r0
	     //now r0 depends on px and py

#ifdef DEBUG
	     std::cout << "r0=" << r0 << " r0inv=" << 1./r0 << " pt=" << pt << std::endl;
#endif
	     //update derivative on D
	     dAPdpx = -(r-r0)*invcurvature*px*pt2inv;//r0 is now 1./r0 (this could go above the redefinition of r0!)
	     dAPdpy = -(r-r0)*invcurvature*py*pt2inv;
	     r0 = 1./r0;//WARNING, now r0 is r0inv (one less temporary)
	     dAPdx = -x*r0*invcurvature;
	     dAPdy = -y*r0*invcurvature;
	     //reduce temporary variables
	     //dxdx = 1 + k*dAPdx*(px*cosAP - py*sinAP);
	     //dydx = k*dAPdx*(py*cosAP + px*sinAP);
	     //dTDdx -= r0*(x*dxdx + y*dydx);
	     dTDdx -= r0*(x*(1 + k*dAPdx*(px*cosAP - py*sinAP)) + y*(k*dAPdx*(py*cosAP + px*sinAP)));
	     //reuse same temporary variables
	     //dxdy = k*dAPdy*(px*cosAP - py*sinAP);
	     //dydy = 1 + k*dAPdy*(py*cosAP + px*sinAP);
	     //dTDdy -= r0*(x*dxdy + y*dydy);
	     dTDdy -= r0*(x*(k*dAPdy*(px*cosAP - py*sinAP)) + y*(1 + k*dAPdy*(py*cosAP + px*sinAP)));
	     //dxdpx = k*(sinAP + px*cosAP*dAPdpx - py*sinAP*dAPdpx);
	     //dydpx = k*(py*cosAP*dAPdpx + 1. - cosAP + px*sinAP*dAPdpx);
	     //dTDdpx -= r0*(x*dxdpx + y*dydpx);
	     dTDdpx -= r0*(x*(k*(sinAP + px*cosAP*dAPdpx - py*sinAP*dAPdpx)) + y*(k*(py*cosAP*dAPdpx + 1. - cosAP + px*sinAP*dAPdpx)));
	     //dxdpy = k*(px*cosAP*dAPdpy - 1. + cosAP - py*sinAP*dAPdpy);
	     //dydpy = k*(sinAP + py*cosAP*dAPdpy + px*sinAP*dAPdpy);
	     //dTDdpy -= r0*(x*dxdpy + y*(dydpy);
	     dTDdpy -= r0*(x*(k*(px*cosAP*dAPdpy - 1. + cosAP - py*sinAP*dAPdpy)) + y*(k*(sinAP + py*cosAP*dAPdpy + px*sinAP*dAPdpy)));

	  }

	  //std::cout << "dTPdx=" << dTDdx*invcurvature << " dTPdy=" << dTDdy*invcurvature << std::endl;

#ifdef DEBUG
	  std::cout << "iteration end, dump parameters" << std::endl;
	  std::cout << "pos = " << outPar.At(n, 0, 0) << " " << outPar.At(n, 1, 0) << " " << outPar.At(n, 2, 0) << std::endl;
	  std::cout << "mom = " << outPar.At(n, 3, 0) << " " << outPar.At(n, 4, 0) << " " << outPar.At(n, 5, 0) << std::endl;
	  std::cout << "r=" << sqrt( outPar.At(n, 0, 0)*outPar.At(n, 0, 0) + outPar.At(n, 1, 0)*outPar.At(n, 1, 0) ) << " pT=" << sqrt( outPar.At(n, 3, 0)*outPar.At(n, 3, 0) + outPar.At(n, 4, 0)*outPar.At(n, 4, 0) ) << std::endl;
#endif
	}
      
      float& TD=totalDistance;
      float  TP=TD*invcurvature;//totalAngPath
      
#ifdef DEBUG
      std::cout << "TD=" << TD << " TP=" << TP << " arrived at r=" << sqrt(outPar.At(n, 0, 0)*outPar.At(n, 0, 0)+outPar.At(n, 1, 0)*outPar.At(n, 1, 0)) << std::endl;
      std::cout << "pos = " << outPar.At(n, 0, 0) << " " << outPar.At(n, 1, 0) << " " << outPar.At(n, 2, 0) << std::endl;
      std::cout << "mom = " << outPar.At(n, 3, 0) << " " << outPar.At(n, 4, 0) << " " << outPar.At(n, 5, 0) << std::endl;
      float rout = sqrt(outPar.At(n, 0, 0)*outPar.At(n, 0, 0)+outPar.At(n, 1, 0)*outPar.At(n, 1, 0));
      float rin = sqrt(xin*xin+yin*yin);
      float dr2 = (outPar.At(n, 0, 0)-xin)*(outPar.At(n, 0, 0)-xin) + (outPar.At(n, 1, 0)-yin)*(outPar.At(n, 1, 0)-yin);
      std::cout << "cos(TP)=" << cos(TP) << " 1-cosTP=" << 1.-cos(TP) << " dR^2/(2k^2pt^2)=" << (rout-rin)*(rout-rin)/(2.*k*k*pt2) << " dR^2/(2k^2pt^2)_v2=" << dr2/(2.*k*k*pt2) << std::endl;
      std::cout << "sinTP=" << sin(TP) << " dR/(kpt) * sqrt(1-dR^2/(4k^2pt^2))=" << (rout-rin)/(k*pt) * sqrt(1. - (rout-rin)*(rout-rin)/(4.*k*k*pt2)) << std::endl;
#endif

      float& iC=invcurvature;
      float dCdpx = k*pxin*ptinv;
      float dCdpy = k*pyin*ptinv;
      float dTPdx = dTDdx*iC;
      float dTPdy = dTDdy*iC;
      float dTPdpx = (dTDdpx - TD*dCdpx*iC)*iC; // MT change: avoid division
      float dTPdpy = (dTDdpy - TD*dCdpy*iC)*iC; // MT change: avoid division

      //std::cout << "dTPdx=" << dTPdx << " dTPdy=" << dTPdy << " dTPdpx=" << dTPdpx << " dTPdpy=" << dTPdpy << std::endl;
      
      float cosTP, sinTP;
      if (Config::useTrigApprox) {
	sincos4(TP, sinTP, cosTP);
      } else {
	cosTP = cos(TP);
	sinTP = sin(TP);
      }

#ifdef DEBUG
      std::cout 
	<< "sinTP=" << sinTP
	<< " cosTP=" << cosTP
	<< " TD=" << TD
	<< std::endl;
#endif
      //now try to make full jacobian
      //derive these to compute jacobian
      //x = xin + k*(pxin*sinTP-pyin*(1-cosTP));
      //y = yin + k*(pyin*sinTP+pxin*(1-cosTP));
      //z = zin + k*TP*pzin;
      //px = pxin*cosTP-pyin*sinTP;
      //py = pyin*cosTP+pxin*sinTP;
      //pz = pzin;
      //jacobian
      
      errorProp(n,0,0) = 1 + k*dTPdx*(pxin*cosTP - pyin*sinTP);	//dxdx;
      errorProp(n,0,1) = k*dTPdy*(pxin*cosTP - pyin*sinTP);	//dxdy;
      errorProp(n,0,2) = 0.;
      errorProp(n,0,3) = k*(sinTP + pxin*cosTP*dTPdpx - pyin*sinTP*dTPdpx); //dxdpx;
      errorProp(n,0,4) = k*(pxin*cosTP*dTPdpy - 1. + cosTP - pyin*sinTP*dTPdpy);//dxdpy;
      errorProp(n,0,5) = 0.;
      
      errorProp(n,1,0) = k*dTPdx*(pyin*cosTP + pxin*sinTP);	//dydx;
      errorProp(n,1,1) = 1 + k*dTPdy*(pyin*cosTP + pxin*sinTP);	//dydy;
      errorProp(n,1,2) = 0.;
      errorProp(n,1,3) = k*(pyin*cosTP*dTPdpx + 1. - cosTP + pxin*sinTP*dTPdpx);//dydpx;
      errorProp(n,1,4) = k*(sinTP + pyin*cosTP*dTPdpy + pxin*sinTP*dTPdpy); //dydpy;
      errorProp(n,1,5) = 0.;
      
      errorProp(n,2,0) = k*pzin*dTPdx;	//dzdx;
      errorProp(n,2,1) = k*pzin*dTPdy;	//dzdy;
      errorProp(n,2,2) = 1.;
      errorProp(n,2,3) = k*pzin*dTPdpx;//dzdpx;
      errorProp(n,2,4) = k*pzin*dTPdpy;//dzdpy;
      errorProp(n,2,5) = k*TP; //dzdpz;
      
      errorProp(n,3,0) = -dTPdx*(pxin*sinTP + pyin*cosTP);	//dpxdx;
      errorProp(n,3,1) = -dTPdy*(pxin*sinTP + pyin*cosTP);	//dpxdy;
      errorProp(n,3,2) = 0.;
      errorProp(n,3,3) = cosTP - dTPdpx*(pxin*sinTP + pyin*cosTP); //dpxdpx;
      errorProp(n,3,4) = -sinTP - dTPdpy*(pxin*sinTP + pyin*cosTP);//dpxdpy;
      errorProp(n,3,5) = 0.;
      
      errorProp(n,4,0) = -dTPdx*(pyin*sinTP - pxin*cosTP); //dpydx;
      errorProp(n,4,1) = -dTPdy*(pyin*sinTP - pxin*cosTP);	//dpydy;
      errorProp(n,4,2) = 0.;
      errorProp(n,4,3) = +sinTP - dTPdpx*(pyin*sinTP - pxin*cosTP);//dpydpx;
      errorProp(n,4,4) = +cosTP - dTPdpy*(pyin*sinTP - pxin*cosTP);//dpydpy;
      errorProp(n,4,5) = 0.;
      
      errorProp(n,5,0) = 0.;
      errorProp(n,5,1) = 0.;
      errorProp(n,5,2) = 0.;
      errorProp(n,5,3) = 0.;
      errorProp(n,5,4) = 0.;
      errorProp(n,5,5) = 1.;

#ifdef DEBUG
      std::cout << "jacobian iterative" << std::endl;
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,0,0),errorProp(n,0,1),errorProp(n,0,2),errorProp(n,0,3),errorProp(n,0,4),errorProp(n,0,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,1,0),errorProp(n,1,1),errorProp(n,1,2),errorProp(n,1,3),errorProp(n,1,4),errorProp(n,1,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,2,0),errorProp(n,2,1),errorProp(n,2,2),errorProp(n,2,3),errorProp(n,2,4),errorProp(n,2,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,3,0),errorProp(n,3,1),errorProp(n,3,2),errorProp(n,3,3),errorProp(n,3,4),errorProp(n,3,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,4,0),errorProp(n,4,1),errorProp(n,4,2),errorProp(n,4,3),errorProp(n,4,4),errorProp(n,4,5));
      printf("%5f %5f %5f %5f %5f %5f\n", errorProp(n,5,0),errorProp(n,5,1),errorProp(n,5,2),errorProp(n,5,3),errorProp(n,5,4),errorProp(n,5,5));
#endif
    }
}

void applyMaterialEffects(const MPlexQF &hitsRl, const MPlexQF& hitsXi, MPlexLS &outErr, MPlexLV& outPar) {

#pragma simd
  for (int n = 0; n < NN; ++n)
    {
      float radL = hitsRl.ConstAt(n,0,0);
      if (radL<0.0000000000001) continue;//ugly, please fixme
      const float& x = outPar.ConstAt(n,0,0);
      const float& y = outPar.ConstAt(n,0,1);
      const float& px = outPar.ConstAt(n,0,3);//FIXME FOR POLCOORD
      const float& py = outPar.ConstAt(n,0,4);
      const float& pz = outPar.ConstAt(n,0,5);
      float r = sqrt(x*x+y*y);
      float pt = px*px + py*py;
      float p2 = pt + pz*pz;
      pt =sqrt(pt);
      float p = sqrt(p2);
      constexpr float mpi = 0.140; // m=140 MeV, pion
      constexpr float mpi2 = 0.140*0.140; // m=140 MeV, pion
      float beta2 = p2/(p2+mpi2);
      float beta = sqrt(beta2);
      //radiation lenght, corrected for the crossing angle (cos alpha from dot product of radius vector and momentum)
      float invCos = (p*r)/fabs(x*px+y*py);
      radL = radL * invCos; //fixme works only for barrel geom
      // multiple scattering
      // in a reference frame defined by the orthogonal unit vectors: u=(px/p,py/p,pz/p) v=(-py/pt,px/pt,0) s=(-pzpx/pt/p,-pzpy/pt/p,pt/p)
      // we consider two planar angles theta1 and theta2 in the uv and us planes respectively
      // note theta1 and theta2 are different angles but with the same rms value thetaMSC
      // first order approximation: sin_thetaMSC ~ thetaMSC
      // px' = px + (py*p*theta1 + pz*px*theta2)/pt; 
      // py' = py - (px*p*theta1 - pz*py*theta2)/pt;
      // pz' = pz + pt*theta2;
      // this actually changes |p| so that p'^2 = p^2(1+2thetaMSC^2) so we should renormalize everything but we neglect this effect here (we are just inflating uncertainties a bit)
      float thetaMSC = 0.0136*sqrt(radL)*(1.+0.038*log(radL))/(beta*p);// eq 32.15
      float thetaMSC2 = thetaMSC*thetaMSC;
      float thetaMSC2overPt2 = thetaMSC2/(pt*pt);
      outErr.At(n, 3, 3) += (py*py*p*p + pz*pz*px*px)*thetaMSC2overPt2;
      outErr.At(n, 4, 4) += (px*px*p*p + pz*pz*py*py)*thetaMSC2overPt2;
      outErr.At(n, 5, 5) += pt*pt*thetaMSC2;
      outErr.At(n, 3, 4) += -px*py*thetaMSC2;
      outErr.At(n, 3, 5) += -pz*px*thetaMSC2;
      outErr.At(n, 4, 5) += -pz*py*thetaMSC2;
      // std::cout << "beta=" << beta << " p=" << p << std::endl;
      // std::cout << "multiple scattering thetaMSC=" << thetaMSC << " thetaMSC2=" << thetaMSC2 << " radL=" << radL << " cxx=" << (py*py*p*p + pz*pz*px*px)*thetaMSC2overPt2 << " cyy=" << (px*px*p*p + pz*pz*py*py)*thetaMSC2overPt2 << " czz=" << pt*pt*thetaMSC2 << std::endl;
      // energy loss
      float gamma = 1./sqrt(1 - beta2);
      float gamma2 = gamma*gamma;
      constexpr float me = 0.0005; // m=0.5 MeV, electron
      float wmax = 2.*me*beta2*gamma2 / ( 1 + 2.*gamma*me/mpi + me*me/(mpi*mpi) );
      constexpr float I = 16.0e-9 * 10.75;
      float deltahalf = log(28.816e-9 * sqrt(2.33*0.498)/I) + log(beta*gamma) - 0.5;
      float dEdx = hitsXi.ConstAt(n,0,0) * invCos * (0.5*log(2*me*beta2*gamma2*wmax/(I*I)) - beta2 - deltahalf) / beta2 ;
      dEdx = dEdx*2.;//xi in cmssw is defined with an extra factor 0.5 with respect to formula 27.1 in pdg
      // std::cout << "dEdx=" << dEdx << " delta=" << deltahalf << std::endl;
      float dP = dEdx/beta;
      outPar.At(n, 0, 3) -= dP*px/p;
      outPar.At(n, 0, 4) -= dP*py/p;
      outPar.At(n, 0, 5) -= dP*pz/p;
      //assume 100% uncertainty
      dP = dP*dP;//warning, redefining dP!
      p2 = 1./p2;//warning, redefining p2!
      outErr.At(n, 3, 3) += dP*px*px*p2;//dP^2*px*px/p^2
      outErr.At(n, 4, 4) += dP*py*py*p2;
      outErr.At(n, 5, 5) += dP*pz*pz*p2;
      p2 = p2/p;//warning, redefining p2!
      outErr.At(n, 3, 4) += dP*px*py*p2;//dP^2*px*py/p^3
      outErr.At(n, 3, 5) += dP*px*pz*p2;
      outErr.At(n, 4, 5) += dP*pz*py*p2;
    }

}

void propagateHelixToRMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const MPlexHV& msPar, 
			          MPlexLS &outErr,       MPlexLV& outPar)
{
#ifdef DEBUG
  const bool dump = false;
#endif

   const idx_t N  = NN;

   outErr = inErr;
   outPar = inPar;

   MPlexLL errorProp;

   MPlexQF msRad;
   // MPlexQF hitsRl;
   // MPlexQF hitsXi;
#pragma simd
   for (int n = 0; n < N; ++n) {
     msRad.At(n, 0, 0) = hipo(msPar.ConstAt(n, 0, 0), msPar.ConstAt(n, 1, 0));
     // if (Config::useCMSGeom) {
     //   hitsRl.At(n, 0, 0) = getRlVal(msRad.ConstAt(n, 0, 0), outPar.ConstAt(n, 2, 0));
     //   hitsXi.At(n, 0, 0) = getXiVal(msRad.ConstAt(n, 0, 0), outPar.ConstAt(n, 2, 0));
     // }
   }

#ifdef POLCOORD
   helixAtRFromIterativePolar(inPar, inChg, outPar, msRad, errorProp);
#else
   helixAtRFromIterative(inPar, inChg, outPar, msRad, errorProp);
#endif

#ifdef DEBUG
   if (dump) {
     for (int kk = 0; kk < N; ++kk)
     {
       printf("outErr before prop %d\n", kk);
       for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
           printf("%8f ", outErr.At(kk,i,j)); printf("\n");
       } printf("\n");

       printf("errorProp %d\n", kk);
       for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
           printf("%8f ", errorProp.At(kk,i,j)); printf("\n");
       } printf("\n");

     }
   }
#endif

   // Matriplex version of:
   // result.errors = ROOT::Math::Similarity(errorProp, outErr);
   MPlexLL temp;
   MultHelixProp      (errorProp, outErr, temp);
   MultHelixPropTransp(errorProp, temp,   outErr);

   // if (Config::useCMSGeom) {
   //   applyMaterialEffects(hitsRl, hitsXi, outErr, outPar);
   // }

   // This dump is now out of its place as similarity is done with matriplex ops.
#ifdef DEBUG
   if (dump) {
     for (int kk = 0; kk < N; ++kk)
     {
       printf("outErr %d\n", kk);
       for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
           printf("%8f ", outErr.At(kk,i,j)); printf("\n");
       } printf("\n");

       printf("outPar %d\n", kk);
       for (int i = 0; i < 6; ++i) {
           printf("%8f ", outPar.At(kk,i,0)); printf("\n");
       } printf("\n");
     }
   }
#endif

#ifdef DEBUG
   if (fabs(hipo(outPar.At(0,0,0), outPar.At(0,1,0))-hipo(msPar.ConstAt(0, 0, 0), msPar.ConstAt(0, 1, 0)))>0.0001) {
     std::cout << "DID NOT GET TO R, dR=" << fabs(hipo(outPar.At(0,0,0), outPar.At(0,1,0))-hipo(msPar.ConstAt(0, 0, 0), msPar.ConstAt(0, 1, 0)))
	       << " r=" << hipo(msPar.ConstAt(0, 0, 0), msPar.ConstAt(0, 1, 0)) << " r0in=" << hipo(inPar.ConstAt(0,0,0), inPar.ConstAt(0,1,0)) << " rout=" << hipo(outPar.At(0,0,0), outPar.At(0,1,0)) << std::endl;
     std::cout << "pt=" << hipo(inPar.ConstAt(0,3,0), inPar.ConstAt(0,4,0)) << " pz=" << inPar.ConstAt(0,5,0) << std::endl;
   }
#endif
}

void propagateHelixToRMPlex(const MPlexLS& inErr,  const MPlexLV& inPar,
                            const MPlexQI& inChg,  const float    r,
			    MPlexLS&       outErr, MPlexLV&       outPar,
                            const int      N_proc)
{
#ifdef DEBUG
  const bool dump = false;
#endif

   outErr = inErr;
   outPar = inPar;

   MPlexLL errorProp;


   MPlexQF msRad;
#pragma simd
   for (int n = 0; n < N_proc; ++n) {
     msRad.At(n, 0, 0) = r;
   }

#ifdef POLCOORD
   helixAtRFromIterativePolar(inPar, inChg, outPar, msRad, errorProp);
#else
   helixAtRFromIterative(inPar, inChg, outPar, msRad, errorProp);
#endif

   //add multiple scattering uncertainty and energy loss (FIXME: in this way it is not applied in track fit)
   if (Config::useCMSGeom) {
     MPlexQF hitsRl;
     MPlexQF hitsXi;
#pragma simd
     for (int n = 0; n < N_proc; ++n) {
       hitsRl.At(n, 0, 0) = getRlVal(r, outPar.ConstAt(n, 2, 0));
       hitsXi.At(n, 0, 0) = getXiVal(r, outPar.ConstAt(n, 2, 0));
     }
     applyMaterialEffects(hitsRl, hitsXi, outErr, outPar);
   }

   // Matriplex version of:
   // result.errors = ROOT::Math::Similarity(errorProp, outErr);

   //MultHelixProp can be optimized for polar coordinates, see GenMPlexOps.pl
   MPlexLL temp;
   MultHelixProp      (errorProp, outErr, temp);
   MultHelixPropTransp(errorProp, temp,   outErr);
   
   // This dump is now out of its place as similarity is done with matriplex ops.
#ifdef DEBUG
   if (dump) {
     for (int kk = 0; kk < N_proc; ++kk)
     {
       printf("outErr %d\n", kk);
       for (int i = 0; i < 6; ++i) { for (int j = 0; j < 6; ++j)
           printf("%8f ", outErr.At(kk,i,j)); printf("\n");
       } printf("\n");

       printf("outPar %d\n", kk);
       for (int i = 0; i < 6; ++i) {
           printf("%8f ", outPar.At(kk,i,0)); printf("\n");
       } printf("\n");
     }
   }
#endif

   /*
     if (fabs(sqrt(outPar[0]*outPar[0]+outPar[1]*outPar[1])-r)>0.0001) {
     std::cout << "DID NOT GET TO R, dR=" << fabs(sqrt(outPar[0]*outPar[0]+outPar[1]*outPar[1])-r)
     << " r=" << r << " r0in=" << r0in << " rout=" << sqrt(outPar[0]*outPar[0]+outPar[1]*outPar[1]) << std::endl;
     std::cout << "pt=" << pt << " pz=" << inPar.At(n, 2) << std::endl;
     }
   */
}
