#include "ConformalUtilsMPlex.h"

inline
void CFMap(const MPlexHH& A, const MPlexHV& B, MPlexHV& C)
{
  using idx_t = Matriplex::idx_t;

  // C = A * B, C is 3x1, A is 3x3 , B is 3x1

  typedef float T;
  typedef float Tv;
  const idx_t N = NN;

  const T *a = A.fArray; ASSUME_ALIGNED(a, 64);
  const Tv *b = B.fArray; ASSUME_ALIGNED(b, 64);
  Tv *c = C.fArray; ASSUME_ALIGNED(c, 64);

#include "CFMatrix33Vector3.ah"
}

//M. Hansroul, H. Jeremie and D. Savard, NIM A 270 (1988) 498
//http://www.sciencedirect.com/science/article/pii/016890028890722X

void conformalFitMPlex(bool fitting, const MPlexQI inChg, 
		       MPlexLS& outErr, MPlexLV& outPar, 
		       const MPlexHV& msPar0, const MPlexHV& msPar1, const MPlexHV& msPar2)
{
  using idx_t = Matriplex::idx_t;
  const idx_t N = NN;

  // Store positions in mplex vectors... could consider storing in a 3x3 matrix, too
  MPlexHV x, y, z, r2;
#pragma simd
  for (int n = 0; n < N; ++n) 
  {
    x.At(n, 0, 0) = msPar0.ConstAt(n, 0, 0);
    x.At(n, 1, 0) = msPar1.ConstAt(n, 0, 0);
    x.At(n, 2, 0) = msPar2.ConstAt(n, 0, 0);

    y.At(n, 0, 0) = msPar0.ConstAt(n, 1, 0);
    y.At(n, 1, 0) = msPar1.ConstAt(n, 1, 0);
    y.At(n, 2, 0) = msPar2.ConstAt(n, 1, 0);

    z.At(n, 0, 0) = msPar0.ConstAt(n, 2, 0);
    z.At(n, 1, 0) = msPar1.ConstAt(n, 2, 0);
    z.At(n, 2, 0) = msPar2.ConstAt(n, 2, 0);
    
    for (int i = 0; i < 3; ++i)
    {
      r2.At(n, i, 0) = getRad2(x.ConstAt(n, i, 0), y.ConstAt(n, i, 0));
    }
  }
  
  MPlexQF initPhi;
  MPlexQI xtou; // bool to determine "split space", i.e. map x to u or v
#pragma simd
  for (int n = 0; n < N; ++n) 
  {
    initPhi.At(n, 0, 0) = std::abs(getPhi(x.ConstAt(n, 0, 0), y.ConstAt(n, 0, 0)));
    xtou.At(n, 0, 0)    = ((initPhi.ConstAt(n, 0, 0) < Config::PIOver4 || initPhi.ConstAt(n, 0, 0) > Config::PI3Over4) ? 1 : 0);
  }

  MPlexHV u,v;
#pragma simd
  for (int n = 0; n < N; ++n) 
  {
    if (xtou.At(n, 0, 0)) // x mapped to u
    {
      for (int i = 0; i < 3; ++i) 
      {
	u.At(n, i, 0) = x.ConstAt(n, i, 0) / r2.ConstAt(n, i, 0);
	v.At(n, i, 0) = y.ConstAt(n, i, 0) / r2.ConstAt(n, i, 0);
      }
    }
    else // x mapped to v
    {
      for (int i = 0; i < 3; ++i) 
      {
	u.At(n, i, 0) = y.ConstAt(n, i, 0) / r2.ConstAt(n, i, 0);
	v.At(n, i, 0) = x.ConstAt(n, i, 0) / r2.ConstAt(n, i, 0);
      }
    }
  }

  MPlexHH A;
  MPlexHV B;
#pragma simd
  for (int n = 0; n < N; ++n) 
  {
    for (int i = 0; i < 3; ++i) 
    {
      A.At(n, i, 0) = 1.0f;
      A.At(n, i, 1) = -u.ConstAt(n, i, 0);
      A.At(n, i, 2) = -u.ConstAt(n, i, 0)*u.ConstAt(n, i, 0);
      B.At(n, i, 0) = v.ConstAt(n, i, 0);
    }
  }
  Matriplex::InvertCramer(A);  
  MPlexHV C; 
  CFMap(A, B, C);
  
  MPlexQF a,b;
#pragma simd
  for (int n = 0; n < N; ++n) 
  {
    b.At(n, 0, 0) = 1.0f/(2.0f*C.ConstAt(n, 0, 0)); 
    a.At(n, 0, 0) = b.ConstAt(n, 0, 0)*C.ConstAt(n, 1, 0); 
  }  

  // do i really need all these temp mplexs????
  MPlexQF pT, pT2, px, py, pz, pz2;
  //#pragma simd
  for (int n = 0; n < N; ++n)
  {
    const float vrx = (xtou.ConstAt(n, 0, 0) ? x.ConstAt(n, 0, 0) - a.ConstAt(n, 0, 0) : x.ConstAt(n, 0, 0) - b.ConstAt(n, 0, 0));
    const float vry = (xtou.ConstAt(n, 0, 0) ? y.ConstAt(n, 0, 0) - b.ConstAt(n, 0, 0) : y.ConstAt(n, 0, 0) - a.ConstAt(n, 0, 0));
    const float invvrxy = 1.0f/std::sqrt(vrx*vrx + vry*vry);
    const float cosphi = vry*invvrxy;
    const float sinphi = vrx*invvrxy;
    const float pT = (-Config::sol*Config::Bfield)*hipo(vrx, vry) / (inChg.ConstAt(n, 0, 0) * 100);
    px.At (n, 0, 0) = std::copysign(pT * cosphi, x.ConstAt(n, 1, 0) - x.ConstAt(n, 0, 0));
    py.At (n, 0, 0) = std::copysign(pT * sinphi, y.ConstAt(n, 1, 0) - y.ConstAt(n, 0, 0));
    const float tantheta = hipo(x.ConstAt(n, 2, 0) - x.ConstAt(n, 0, 0), y.ConstAt(n, 2, 0) - y.ConstAt(n, 0, 0))/(z.ConstAt(n, 2, 0) - z.ConstAt(n, 0, 0));
    pz.At (n, 0, 0) = std::copysign(pT / tantheta, z.ConstAt(n, 1, 0) - z.ConstAt(n, 0, 0));

    pT2.At(n, 0, 0) = pT*pT;
    pz2.At(n, 0, 0) = pz.ConstAt(n, 0, 0)*pz.ConstAt(n, 0, 0);
  }

#ifdef INWARDFIT
  if (fitting)
  {
#pragma simd
    for (int n = 0; n < N; ++n)
    {
      px.At(n, 0, 0) *= -1.0f;
      py.At(n, 0, 0) *= -1.0f;
      pz.At(n, 0, 0) *= -1.0f;
    }
  }
#endif
 
 //  Start setting the output parameters
#pragma simd
  for (int n = 0; n < N; ++n)
  {
    outPar.At(n, 0, 0) = x.ConstAt(n, 0, 0);
    outPar.At(n, 1, 0) = y.ConstAt(n, 0, 0);
    outPar.At(n, 2, 0) = z.ConstAt(n, 0, 0);
    outPar.At(n, 3, 0) = px.ConstAt(n, 0, 0);
    outPar.At(n, 4, 0) = py.ConstAt(n, 0, 0);
    outPar.At(n, 5, 0) = pz.ConstAt(n, 0, 0);
  }

  // Use r-phi smearing to set initial error estimation
  // uncertainties set by hand by making pulls width of 1.0 (extracted from residuals)
  float ptinverr, phierr, thetaerr;
  if (fitting)
  {
    ptinverr = Config::ptinverr049;
    phierr   = Config::phierr049;
    thetaerr = Config::thetaerr049;
  }
  else
  {
    ptinverr = Config::ptinverr012;
    phierr   = Config::phierr012;
    thetaerr = Config::thetaerr012;
  }

  MPlexQF varPt, varPhi, invvarR2, varTheta;
#pragma simd
  for (int n = 0; n < N; ++n)
  {
    varPt.At   (n, 0, 0) = pT2.ConstAt(n, 0, 0)*ptinverr*ptinverr;
    varPhi.At  (n, 0, 0) = Config::varXY/r2.ConstAt(n, 0, 0);
    invvarR2.At(n, 0, 0) = Config::varR/r2.ConstAt(n, 0, 0);
    varTheta.At(n, 0, 0) = ((pT2.ConstAt(n, 0, 0) + pz2.ConstAt(n, 0, 0))*(pT2.ConstAt(n, 0, 0) + pz2.ConstAt(n, 0, 0))) / pT2.ConstAt(n, 0, 0) * thetaerr * thetaerr;
  }

#pragma simd
  for (int n = 0; n < N; ++n)
  {
    outErr.At(n, 0, 0) = x.ConstAt(n, 0, 0)*x.ConstAt(n, 0, 0)*invvarR2.ConstAt(n, 0, 0) + y.ConstAt(n, 0, 0)*y.ConstAt(n, 0, 0)*varPhi.ConstAt(n, 0, 0);
    outErr.At(n, 0, 1) = x.ConstAt(n, 0, 0)*y.ConstAt(n, 0, 0)*(invvarR2.ConstAt(n, 0, 0) - varPhi.ConstAt(n, 0, 0));
    outErr.At(n, 0, 2) = 0.;
    outErr.At(n, 0, 3) = 0.;
    outErr.At(n, 0, 4) = 0.;
    outErr.At(n, 0, 5) = 0.;

    outErr.At(n, 1, 0) = outErr.ConstAt(n, 0, 1);
    outErr.At(n, 1, 1) = y.ConstAt(n, 0, 0)*y.ConstAt(n, 0, 0)*invvarR2.ConstAt(n, 0, 0) + x.ConstAt(n, 0, 0)*x.ConstAt(n, 0, 0)*varPhi.ConstAt(n, 0, 0);
    outErr.At(n, 1, 2) = 0.;
    outErr.At(n, 1, 3) = 0.;
    outErr.At(n, 1, 4) = 0.;
    outErr.At(n, 1, 5) = 0.;

    outErr.At(n, 2, 0) = 0.;
    outErr.At(n, 2, 1) = 0.;
    outErr.At(n, 2, 2) = Config::varZ;
    outErr.At(n, 2, 3) = 0.;
    outErr.At(n, 2, 4) = 0.;
    outErr.At(n, 2, 5) = 0.;

    outErr.At(n, 3, 0) = 0.;
    outErr.At(n, 3, 1) = 0.;
    outErr.At(n, 3, 2) = 0.;
    outErr.At(n, 3, 3) = px.ConstAt(n, 0, 0)*px.ConstAt(n, 0, 0)*varPt(n, 0, 0) + py.ConstAt(n, 0, 0)*py.ConstAt(n, 0, 0)*phierr*phierr;
    outErr.At(n, 3, 4) = 0.;
    outErr.At(n, 3, 5) = 0.;

    outErr.At(n, 4, 0) = 0.;
    outErr.At(n, 4, 1) = 0.;
    outErr.At(n, 4, 2) = 0.;
    outErr.At(n, 4, 3) = 0.;
    outErr.At(n, 4, 4) = py.ConstAt(n, 0, 0)*py.ConstAt(n, 0, 0)*varPt(n, 0, 0) + px.ConstAt(n, 0, 0)*px.ConstAt(n, 0, 0)*phierr*phierr;
    outErr.At(n, 4, 5) = 0.;

    outErr.At(n, 5, 0) = 0.;
    outErr.At(n, 5, 1) = 0.;
    outErr.At(n, 5, 2) = 0.;
    outErr.At(n, 5, 3) = 0.;
    outErr.At(n, 5, 4) = 0.;
    outErr.At(n, 5, 5) = pz2.ConstAt(n, 0, 0)*varPt.ConstAt(n, 0, 0) + varTheta.ConstAt(n, 0, 0);
  }  
}
