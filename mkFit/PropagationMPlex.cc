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

// this version does not assume to know which elements are 0 or 1, so it does the full mupltiplication
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

}



void computeJacobianSimple(int n, MPlexLL& errorProp, float s, float k, float p, float pxin, float pyin, float pzin, float TP, float cosTP, float sinTP) {

  // std::cout << "total path s=" << s << std::endl;
  // TD = s*pt/p;
  // TP = TD/(pt*k) = s/(p*k);
  float dTPdpx = -s*pxin/(k*p*p*p);
  float dTPdpy = -s*pyin/(k*p*p*p);
  float dTPdpz = -s*pzin/(k*p*p*p);
  //ok let's assume that the quantity with no error is the angular path (phase change)
  dTPdpx = 0;
  dTPdpy = 0;
  dTPdpz = 0;
  
  //derive these to compute jacobian
  //x = xin + k*(pxin*sinTP-pyin*(1-cosTP));
  //y = yin + k*(pyin*sinTP+pxin*(1-cosTP));
  //z = zin + k*TP*pzin;
  //px = pxin*cosTP-pyin*sinTP;
  //py = pyin*cosTP+pxin*sinTP;
  //pz = pzin;
  //jacobian
  
  errorProp(n,0,0) = 1.;	                                             //dxdx
  errorProp(n,0,1) = 0.;	                                             //dxdy
  errorProp(n,0,2) = 0.;                                                     //dxdz
  errorProp(n,0,3) = k*(sinTP + pxin*cosTP*dTPdpx - pyin*sinTP*dTPdpx);      //dxdpx
  errorProp(n,0,4) = k*(pxin*cosTP*dTPdpy - 1. + cosTP - pyin*sinTP*dTPdpy); //dxdpy
  errorProp(n,0,5) = k*dTPdpz*(pxin*cosTP - pyin*sinTP);                     //dxdpz
  errorProp(n,1,0) = 0.;	                                             //dydx
  errorProp(n,1,1) = 1.;	                                             //dydy
  errorProp(n,1,2) = 0.;                                                     //dydz
  errorProp(n,1,3) = k*(pyin*cosTP*dTPdpx + 1. - cosTP + pxin*sinTP*dTPdpx); //dydpx
  errorProp(n,1,4) = k*(sinTP + pyin*cosTP*dTPdpy + pxin*sinTP*dTPdpy);      //dydpy
  errorProp(n,1,5) = k*dTPdpz*(pyin*cosTP + pxin*sinTP);                     //dydpz
  errorProp(n,2,0) = 0.;	                                             //dzdx
  errorProp(n,2,1) = 0.;	                                             //dzdy
  errorProp(n,2,2) = 1.;                                                     //dzdz
  errorProp(n,2,3) = k*pzin*dTPdpx;                                          //dzdpx
  errorProp(n,2,4) = k*pzin*dTPdpy;                                          //dzdpy
  errorProp(n,2,5) = k*(TP + dTPdpz*pzin);                                   //dzdpz
  errorProp(n,3,0) = 0.;	                                             //dpxdx
  errorProp(n,3,1) = 0.;	                                             //dpxdy
  errorProp(n,3,2) = 0.;                                                     //dpxdz
  errorProp(n,3,3) = cosTP - dTPdpx*(pxin*sinTP + pyin*cosTP);               //dpxdpx
  errorProp(n,3,4) = -sinTP - dTPdpy*(pxin*sinTP + pyin*cosTP);              //dpxdpy
  errorProp(n,3,5) = -dTPdpz*(pxin*sinTP + pyin*cosTP);                      //dpxdpz
  errorProp(n,4,0) = 0.;                                                     //dpydx
  errorProp(n,4,1) = 0.;	                                             //dpydy
  errorProp(n,4,2) = 0.;                                                     //dpydz
  errorProp(n,4,3) = +sinTP - dTPdpx*(pyin*sinTP - pxin*cosTP);              //dpydpx
  errorProp(n,4,4) = +cosTP - dTPdpy*(pyin*sinTP - pxin*cosTP);              //dpydpy
  errorProp(n,4,5) = -dTPdpz*(pyin*sinTP - pxin*cosTP);                      //dpydpz
  errorProp(n,5,0) = 0.;                                                     //dpzdx
  errorProp(n,5,1) = 0.;						     //dpzdy
  errorProp(n,5,2) = 0.;						     //dpzdz 
  errorProp(n,5,3) = 0.;						     //dpzdpx
  errorProp(n,5,4) = 0.;						     //dpzdpy
  errorProp(n,5,5) = 1.;						     //dpzdpz  
  
}

//#define DEBUG

void helixAtRFromIterative(const MPlexLV& inPar, const MPlexQI& inChg, MPlexLV& outPar, const MPlexQF &msRad, MPlexLL& errorProp) {

#pragma simd
  for (int n = 0; n < NN; ++n)
    {
      const float& xin  = inPar.ConstAt(n, 0, 0);
      const float& yin  = inPar.ConstAt(n, 1, 0);
      const float& pxin = inPar.ConstAt(n, 3, 0);
      const float& pyin = inPar.ConstAt(n, 4, 0);
      const float& pzin = inPar.ConstAt(n, 5, 0);
      const float& r    = msRad.ConstAt(n, 0, 0);
      float r0 = hipo(xin, yin);
      
#ifdef DEBUG
      std::cout << "attempt propagation from r=" << r0 << " to r=" << r << std::endl;
      std::cout << "x=" << xin << " y=" << yin  << " z=" << inPar.ConstAt(n, 2, 0) << " px=" << pxin << " py=" << pyin << " pz=" << pzin << " q=" << inChg.ConstAt(n, 0, 0) << std::endl;
      // if ((r0-r)>=0) {
      //    if (dump) std::cout << "target radius same or smaller than starting point, returning input" << std::endl;
      //    return;
      // }
#endif

      if (fabs(r-r0)<0.0001) {
#ifdef DEBUG
	std::cout << "distance less than 1mum, skip" << std::endl;
#endif
	computeJacobianSimple(n, errorProp, 0, 1, 1, 1, 1, 1, 0, 1, 0);//get an identity matrix
	continue;
      }
      
      float pt2    = pxin*pxin+pyin*pyin;
      float pt     = sqrt(pt2);
      float ptinv  = 1./pt;
      float pt2inv = ptinv*ptinv;
      //p=0.3Br => r=p/(0.3*B)
      float k = inChg.ConstAt(n, 0, 0) * 100. / (-0.299792458*Config::Bfield);
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
      // float dxdvar = 0.;
      // float dydvar = 0.;
      //5 iterations is a good starting point
      //const unsigned int Niter = 10;
      // const unsigned int Niter = 5+std::round(r-r0)/2;
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

	  if (Config::useSimpleJac==0 && i+1 != Config::Niter && r0 > 0 && fabs((r-r0)*invcurvature)>0.000000001)
         {
            //update derivatives on total distance for next step, where totalDistance+=r-r0
            //now r0 depends on px and py
            r0 = 1./r0;//WARNING, now r0 is r0inv (one less temporary)

#ifdef DEBUG
            std::cout << "r0=" << 1./r0 << " r0inv=" << r0 << " pt=" << pt << std::endl;
#endif

            //update derivative on D
            dAPdx = -x*r0*invcurvature;
            dAPdy = -y*r0*invcurvature;
            dAPdpx = -(r-1./r0)*invcurvature*px*pt2inv;//weird, using r0 instead of 1./r0 improves things but it should be wrong since r0 in now r0inv
            dAPdpy = -(r-1./r0)*invcurvature*py*pt2inv;//weird, using r0 instead of 1./r0 improves things but it should be wrong since r0 in now r0inv
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
            //dTDdpy -= r0*(x*dxdpy + y*(k*dydpy);
            dTDdpy -= r0*(x*(k*(px*cosAP*dAPdpy - 1. + cosAP - py*sinAP*dAPdpy)) + y*(k*(sinAP + py*cosAP*dAPdpy + px*sinAP*dAPdpy)));

         }
	  
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
#endif

      float& iC=invcurvature;
      float dCdpx = k*pxin*ptinv;
      float dCdpy = k*pyin*ptinv;
      float dTPdx = dTDdx*iC;
      float dTPdy = dTDdy*iC;
      float dTPdpx = (dTDdpx - TD*dCdpx*iC)*iC; // MT change: avoid division
      float dTPdpy = (dTDdpy - TD*dCdpy*iC)*iC; // MT change: avoid division
      
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

      if (Config::useSimpleJac) { 
	//assume total path length s as given and with no uncertainty
	float p = pt2 + pzin*pzin;
	p = sqrt(p);
	float s = TD*p*ptinv;
	computeJacobianSimple(n, errorProp, s, k, p, pxin, pyin, pzin, TP, cosTP, sinTP);
      } else {
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
      }

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

void helixAtRFromIntersection(const MPlexLV& inPar, const MPlexQI& inChg, MPlexLV& outPar, const MPlexQF &msRad, MPlexLL& errorProp) {

#pragma simd
  for (int n = 0; n < NN; ++n)
    {
      const float& xin  = inPar.ConstAt(n, 0, 0);
      const float& yin  = inPar.ConstAt(n, 1, 0);
      const float& zin  = inPar.ConstAt(0, 2, 0);
      const float& pxin = inPar.ConstAt(n, 3, 0);
      const float& pyin = inPar.ConstAt(n, 4, 0);
      const float& pzin = inPar.ConstAt(n, 5, 0);
      const float& rout = msRad.ConstAt(n, 0, 0);
      const float r0in = hipo(xin, yin);
      
#ifdef DEBUG
      std::cout << "attempt propagation from r=" << r0in << " to r=" << rout << std::endl;
      std::cout << "x=" << xin << " y=" << yin  << " z=" << inPar.ConstAt(n, 2, 0) << " px=" << pxin << " py=" << pyin << " pz=" << pzin << " q=" << inChg.ConstAt(n, 0, 0) << std::endl;
      // if ((r0in-r)>=0) {
      //    if (dump) std::cout << "target radius same or smaller than starting point, returning input" << std::endl;
      //    return;
      // }
#endif

      if (fabs(rout-r0in)<0.0001) {
#ifdef DEBUG
	std::cout << "distance less than 1mum, skip" << std::endl;
#endif
	computeJacobianSimple(n, errorProp, 0, 1, 1, 1, 1, 1, 0, 1, 0);//get an identity matrix
	continue;
      }
      
      float pt2    = pxin*pxin+pyin*pyin;
      float pt     = sqrt(pt2);
      float ptinv  = 1./pt;
      float pt2inv = ptinv*ptinv;
      //p=0.3Br => r=p/(0.3*B)
      float k = inChg.ConstAt(n, 0, 0) * 100. / (-0.299792458*Config::Bfield);
      float curvature = pt*k;//in cm
      float invcurvature = 1./(curvature);//in 1./cm
      
#ifdef DEBUG
      std::cout << "k=" << k << " curvature=" << curvature << " invcurvature=" << invcurvature << std::endl;
#endif
      
      //coordinates of center of circle defined as helix projection on transverse plane 
      float xc = xin - k*pyin;
      float yc = yin + k*pxin;
      float rc = sqrt(xc*xc+yc*yc);
      
      float TP, x, y = 0.;
      if (fabs(xc)>fabs(yc) || fabs(yc)>0.001) {
	//solve for x since yc!=0
	float A = 1 + xc*xc/(yc*yc);
	float B = -(rout*rout + rc*rc - curvature*curvature)*xc/(yc*yc);
	float C = (rout*rout + rc*rc - curvature*curvature)*(rout*rout + rc*rc - curvature*curvature)/(4*yc*yc) - rout*rout;
	//first solution
	float x_p = ( -B + sqrt(B*B - 4*A*C) ) / (2*A);
	float y_p = -x_p*xc/yc + (rout*rout + rc*rc - curvature*curvature)/(2*yc);
	float cosDelta_p = (pxin*(x_p-xin) + pyin*(y_p-yin))/(pt*sqrt((x_p-xin)*(x_p-xin)+(y_p-yin)*(y_p-yin)));
	//second solution
	float x_m = ( -B - sqrt(B*B - 4*A*C) ) / (2*A);
	float y_m = -x_m*xc/yc + (rout*rout + rc*rc - curvature*curvature)/(2*yc);
	float cosDelta_m = (pxin*(x_m-xin) + pyin*(y_m-yin))/(pt*sqrt((x_m-xin)*(x_m-xin)+(y_m-yin)*(y_m-yin)));
#ifdef DEBUG
	std::cout << "solve for x" << std::endl;
	std::cout << "xc=" << xc << " yc=" << yc << " rc=" << rc << std::endl;
	std::cout << "A=" << A << " B=" << B << " C=" << C << std::endl;
	std::cout << "xp=" << x_p << " y_p=" << y_p << " cosDelta_p=" << cosDelta_p << std::endl;
	std::cout << "xm=" << x_m << " y_m=" << y_m << " cosDelta_m=" << cosDelta_m << std::endl;
#endif 
	float Sx = sqrt(B*B - 4*A*C);
	//arbitrate based on momentum and vector connecting the end points
	if ( (rout>r0in) ? (cosDelta_p > cosDelta_m) : (cosDelta_p < cosDelta_m)) { 
	  float chord_p = sqrt( (x_p-xin)*(x_p-xin) + (y_p-yin)*(y_p-yin) );
	  float sinTPHalf_p = 0.5*chord_p*invcurvature;
	  TP = 2*asin(sinTPHalf_p);
	  x = x_p;
	  y = y_p;
#ifdef DEBUG
	  std::cout << "pos solution" << std::endl;
	  std::cout << "chord_p=" << chord_p << " sinTPHalf_p=" << sinTPHalf_p << " TP=" << TP << std::endl;
#endif
	} else {
	  float chord_m = sqrt( (x_m-xin)*(x_m-xin) + (y_m-yin)*(y_m-yin) );
	  float sinTPHalf_m = 0.5*chord_m*invcurvature;
	  TP = 2*asin(sinTPHalf_m);
	  x = x_m;
	  y = y_m;
#ifdef DEBUG
	  std::cout << "neg solution" << std::endl;
	  std::cout << "chord_m=" << chord_m << " sinTPHalf_m=" << sinTPHalf_m << " TP=" << TP << std::endl;
#endif
	} 
      } else {
	//solve for y since xc!=0
	float A = 1 + yc*yc/(xc*xc);
	float B = -(rout*rout + rc*rc - curvature*curvature)*yc/(xc*xc);
	float C = (rout*rout + rc*rc - curvature*curvature)*(rout*rout + rc*rc - curvature*curvature)/(4*xc*xc) - rout*rout;
	//first solution
	float y_p = ( -B + sqrt(B*B - 4*A*C) ) / (2*A);
	float x_p = -y_p*yc/xc + (rout*rout + rc*rc - curvature*curvature)/(2*xc);
	float cosDelta_p = (pxin*(x_p-xin) + pyin*(y_p-yin))/(pt*sqrt((x_p-xin)*(x_p-xin)+(y_p-yin)*(y_p-yin)));
	//second solution
	float y_m = ( -B - sqrt(B*B - 4*A*C) ) / (2*A);
	float x_m = -y_m*yc/xc + (rout*rout + rc*rc - curvature*curvature)/(2*xc);
	float cosDelta_m = (pxin*(x_m-xin) + pyin*(y_m-yin))/(pt*sqrt((x_m-xin)*(x_m-xin)+(y_m-yin)*(y_m-yin)));
#ifdef DEBUG
	std::cout << "solve for y" << std::endl;
	std::cout << "xc=" << xc << " yc=" << yc << " rc=" << rc << std::endl;
	std::cout << "A=" << A << " B=" << B << " C=" << C << std::endl;
	std::cout << "xp=" << x_p << " y_p=" << y_p << " cosDelta_p=" << cosDelta_p << std::endl;
	std::cout << "xm=" << x_m << " y_m=" << y_m << " cosDelta_m=" << cosDelta_m << std::endl;
#endif
	float Sx = sqrt(B*B - 4*A*C);
	//arbitrate based on momentum and vector connecting the end points
	if ( (rout>r0in) ? (cosDelta_p > cosDelta_m) : (cosDelta_p < cosDelta_m)) { 
	  float chord_p = sqrt( (x_p-xin)*(x_p-xin) + (y_p-yin)*(y_p-yin) );
	  float sinTPHalf_p = 0.5*chord_p*invcurvature;
	  TP = 2*asin(sinTPHalf_p);
	  x = x_p;
	  y = y_p;
	} else {
	  float chord_m = sqrt( (x_m-xin)*(x_m-xin) + (y_m-yin)*(y_m-yin) );
	  float sinTPHalf_m = 0.5*chord_m*invcurvature;
	  TP = 2*asin(sinTPHalf_m);
	  x = x_m;
	  y = y_m;
	} 
      }

      float cosTP, sinTP;
      if (Config::useTrigApprox) {
	sincos4(TP, sinTP, cosTP);
      } else {
	cosTP = cos(TP);
	sinTP = sin(TP);
      }

      //correct if we got the wrong sign (fixme, find a better way to do it!)
      if ( fabs((xin + k*(pxin*sinTP-pyin*(1-cosTP)))-x)>fabs((xin + k*(pxin*sin(-TP)-pyin*(1-cos(-TP))))-x) ) {
      	   TP=-TP;
      	   cosTP=cos(TP);
      	   sinTP=sin(TP);
      	}

      //helix propagation formulas
      //http://www.phys.ufl.edu/~avery/fitting/fitting4.pdf
      outPar.At(n, 0, 0) = x; //xin + k*(pxin*sinTP-pyin*(1-cosTP));
      outPar.At(n, 1, 0) = y; //yin + k*(pyin*sinTP+pxin*(1-cosTP));
      outPar.At(n, 2, 0) = zin + k*TP*pzin;
      outPar.At(n, 3, 0) = pxin*cosTP-pyin*sinTP;
      outPar.At(n, 4, 0) = pyin*cosTP+pxin*sinTP;
      //outPar.At(n, 5, 0) = pzin; //take this out as it is redundant

#ifdef DEBUG
      std::cout << "propagation end, dump parameters" << std::endl;
      std::cout << "pos = " << outPar.At(n, 0, 0) << " " << outPar.At(n, 1, 0) << " " << outPar.At(n, 2, 0) << std::endl;
      std::cout << "mom = " << outPar.At(n, 3, 0) << " " << outPar.At(n, 4, 0) << " " << outPar.At(n, 5, 0) << std::endl;
      std::cout << "r=" << sqrt( outPar.At(n, 0, 0)*outPar.At(n, 0, 0) + outPar.At(n, 1, 0)*outPar.At(n, 1, 0) ) << " pT=" << sqrt( outPar.At(n, 3, 0)*outPar.At(n, 3, 0) + outPar.At(n, 4, 0)*outPar.At(n, 4, 0) ) << std::endl;
#endif

      float p = pt2 + pzin*pzin;
      p=sqrt(p);
      float s = TP*curvature*p*ptinv;
      computeJacobianSimple(n, errorProp, s, k, p, pxin, pyin, pzin, TP, cosTP, sinTP);

#ifdef DEBUG
      std::cout << "jacobian intersection" << std::endl;
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
      const float& px = outPar.ConstAt(n,0,3);
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
      // px' = px - (py*p*theta1 + pz*px*theta2)/pt; 
      // py' = py + (px*p*theta1 - pz*py*theta2)/pt;
      // pz' = pz + pt*theta2;
      // this actually changes |p| so that p'^2 = p^2(1+2thetaMSC^2) so we should renormalize everything but we neglect this effect here (we are just inflating uncertainties a bit)
      float thetaMSC = 0.0136*sqrt(radL)*(1.+0.038*log(radL))/(beta*p);// eq 32.15
      float thetaMSC2 = thetaMSC*thetaMSC;
      float thetaMSC2overPt2 = thetaMSC2/(pt*pt);
      outErr.At(n, 3, 3) += (py*py*p*p + pz*pz*px*px)*thetaMSC2overPt2;
      outErr.At(n, 4, 4) += (px*px*p*p + pz*pz*py*py)*thetaMSC2overPt2;
      outErr.At(n, 5, 5) += pt*pt*thetaMSC2;
      outErr.At(n, 3, 4) += px*py*(p*p + pz*pz)*thetaMSC2overPt2;
      outErr.At(n, 3, 5) += -pz*px*thetaMSC2;
      outErr.At(n, 4, 5) += -pz*py*thetaMSC2;
      // std::cout << "beta=" << beta << " p=" << p << std::endl;
      std::cout << "multiple scattering thetaMSC=" << thetaMSC << " thetaMSC2=" << thetaMSC2 << " radL=" << radL << " cxx=" << (py*py*p*p + pz*pz*px*px)*thetaMSC2overPt2 << " cyy=" << (px*px*p*p + pz*pz*py*py)*thetaMSC2overPt2 << " czz=" << pt*pt*thetaMSC2 << std::endl;
      // energy loss
      float gamma = 1./sqrt(1 - beta2);
      float gamma2 = gamma*gamma;
      constexpr float me = 0.0005; // m=0.5 MeV, electron
      float wmax = 2.*me*beta2*gamma2 / ( 1 + 2.*gamma*me/mpi + me*me/(mpi*mpi) );
      constexpr float I = 16.0e-9 * 10.75;
      float deltahalf = log(28.816e-9 * sqrt(2.33*0.498)/I) + log(beta*gamma) - 0.5;
      float dEdx = hitsXi.ConstAt(n,0,0) * invCos * (0.5*log(2*me*beta2*gamma2*wmax/(I*I)) - beta2 - deltahalf) / beta2 ;
      // std::cout << "dEdx=" << dEdx << " delta=" << deltahalf << std::endl;
      float dP = dEdx/beta;
      outPar.At(n, 0, 3) -= dP*px/p;
      outPar.At(n, 0, 4) -= dP*py/p;
      outPar.At(n, 0, 5) -= dP*pz/p;
      //we do nothing on the uncertainty for now
    }

}

void propagateHelixToRMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const MPlexHV& msPar, 
			    //const MPlexQF &hitsRl, const MPlexQF& hitsXi,
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
#pragma simd
   for (int n = 0; n < N; ++n) {
     msRad.At(n, 0, 0) = hipo(msPar.ConstAt(n, 0, 0), msPar.ConstAt(n, 1, 0));
   }

   if (Config::doIterative) {
     helixAtRFromIterative(inPar, inChg, outPar, msRad, errorProp);
   } else {
     // std::cout << "FROM ITERATIVE" << std::endl;
     // helixAtRFromIterative(inPar, inChg, outPar, msRad, errorProp);
     // MPlexLV tmpPar = outPar;
     // std::cout << "FROM INTERSECTION" << std::endl;
     // TODO: not ported for now. Assuming Config::doIterative
     helixAtRFromIntersection(inPar, inChg, outPar, msRad, errorProp);
     // if (fabs(tmpPar.At(0,0,0)-outPar.At(0,0,0))>0.01 ) std::cout << "PROPAGATION PROBLEM" << std::endl;
   }

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

   //add multiple scattering uncertainty and energy loss
   //applyMaterialEffects(hitsRl, hitsXi, outErr, outPar);

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

   if (Config::doIterative) {
     helixAtRFromIterative(inPar, inChg, outPar, msRad, errorProp);
   } else {
     helixAtRFromIntersection(inPar, inChg, outPar, msRad, errorProp);
   }

   // Matriplex version of:
   // result.errors = ROOT::Math::Similarity(errorProp, outErr);

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
