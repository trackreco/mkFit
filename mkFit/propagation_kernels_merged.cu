#include "Config.h"
#include "propagation_kernels.h"
#include <stdio.h>

#define L 6
#define LL 36
#define LS 21

// values from 32 to 512 give good results.
// 32 gives slightly better results (on a K40)
#define BLOCK_SIZE_X 32
#define BLOCK_SIZE_X_2 256
#define BLOCK_SIZE_HELIX_AT 128
#define MAX_BLOCKS_X 65535 // CUDA constraint

__device__ float hipo_merged(float x, float y) {
  return sqrt(x*x + y*y);
}
__device__ void sincos4_merged(float x, float& sin, float& cos) {
   // Had this writen with explicit division by factorial.
   // The *whole* fitting test ran like 2.5% slower on MIC, sigh.
   cos  = 1;
   sin  = x;   x *= x * 0.5f;
   cos -= x;   x *= x * 0.33333333f;
   sin -= x;   x *= x * 0.25f;
   cos += x;
}

__device__ void computeJacobianSimple_merged(int n, 
    float *errorProp, size_t errorProp_stride,
    float s, float k, float p, float pxin, float pyin, float pzin, 
    float TP, float cosTP, float sinTP, int N) {

  /*int i = threadIdx.x;*/
  /*size_t epN = errorProp_stride;*/

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
  
  errorProp[(0*L + 0)] = 1.;	                                             //dxdx
  errorProp[(0*L + 1)] = 0.;	                                             //dxdy
  errorProp[(0*L + 2)] = 0.;                                                     //dxdz
  errorProp[(0*L + 3)] = k*(sinTP + pxin*cosTP*dTPdpx - pyin*sinTP*dTPdpx);      //dxdpx
  errorProp[(0*L + 4)] = k*(pxin*cosTP*dTPdpy - 1. + cosTP - pyin*sinTP*dTPdpy); //dxdpy
  errorProp[(0*L + 5)] = k*dTPdpz*(pxin*cosTP - pyin*sinTP);                     //dxdpz
  errorProp[(1*L + 0)] = 0.;	                                             //dydx
  errorProp[(1*L + 1)] = 1.;	                                             //dydy
  errorProp[(1*L + 2)] = 0.;                                                     //dydz
  errorProp[(1*L + 3)] = k*(pyin*cosTP*dTPdpx + 1. - cosTP + pxin*sinTP*dTPdpx); //dydpx
  errorProp[(1*L + 4)] = k*(sinTP + pyin*cosTP*dTPdpy + pxin*sinTP*dTPdpy);      //dydpy
  errorProp[(1*L + 5)] = k*dTPdpz*(pyin*cosTP + pxin*sinTP);                     //dydpz
  errorProp[(2*L + 0)] = 0.;	                                             //dzdx
  errorProp[(2*L + 1)] = 0.;	                                             //dzdy
  errorProp[(2*L + 2)] = 1.;                                                     //dzdz
  errorProp[(2*L + 3)] = k*pzin*dTPdpx;                                          //dzdpx
  errorProp[(2*L + 4)] = k*pzin*dTPdpy;                                          //dzdpy
  errorProp[(2*L + 5)] = k*(TP + dTPdpz*pzin);                                   //dzdpz
  errorProp[(3*L + 0)] = 0.;	                                             //dpxdx
  errorProp[(3*L + 1)] = 0.;	                                             //dpxdy
  errorProp[(3*L + 2)] = 0.;                                                     //dpxdz
  errorProp[(3*L + 3)] = cosTP - dTPdpx*(pxin*sinTP + pyin*cosTP);               //dpxdpx
  errorProp[(3*L + 4)] = -sinTP - dTPdpy*(pxin*sinTP + pyin*cosTP);              //dpxdpy
  errorProp[(3*L + 5)] = -dTPdpz*(pxin*sinTP + pyin*cosTP);                      //dpxdpz
  errorProp[(4*L + 0)] = 0.;                                                     //dpydx
  errorProp[(4*L + 1)] = 0.;	                                             //dpydy
  errorProp[(4*L + 2)] = 0.;                                                     //dpydz
  errorProp[(4*L + 3)] = +sinTP - dTPdpx*(pyin*sinTP - pxin*cosTP);              //dpydpx
  errorProp[(4*L + 4)] = +cosTP - dTPdpy*(pyin*sinTP - pxin*cosTP);              //dpydpy
  errorProp[(4*L + 5)] = -dTPdpz*(pyin*sinTP - pxin*cosTP);                      //dpydpz
  errorProp[(5*L + 0)] = 0.;                                                     //dpzdx
  errorProp[(5*L + 1)] = 0.;						     //dpzdy
  errorProp[(5*L + 2)] = 0.;						     //dpzdz 
  errorProp[(5*L + 3)] = 0.;						     //dpzdpx
  errorProp[(5*L + 4)] = 0.;						     //dpzdpy
  errorProp[(5*L + 5)] = 1.;						     //dpzdpz  
}

/// Compute MsRad /////////////////////////////////////////////////////////////
// Not passing msRad.stride, as QF == 1 (second dim f msRad)
__device__ void computeMsRad_fn(float *msPar, size_t stride_msPar,
    float* msRad, int N) {
  int n = threadIdx.x + blockIdx.x * blockDim.x;
  if (n < N) {
    *msRad = hipo_merged(msPar[n], msPar[n + stride_msPar]);
  }
}

__device__ 
void helixAtRFromIterative_fn(float *inPar, size_t inPar_stride,
    int *inChg, float *outPar, size_t outPar_stride, float msRad, 
    float *sh_errorProp, size_t errorProp_stride, int N) {

  size_t epN = errorProp_stride;
  size_t opN = outPar_stride;
  size_t ipN = inPar_stride;

  int n = threadIdx.x + blockIdx.x * blockDim.x;
  /*int i = threadIdx.x;*/

  float outPar_reg[5];
  /*__shared__ float sh_errorProp[LL][BLOCK_SIZE_HELIX_AT];*/
  /*float sh_errorProp[LL];*/

  if (n < N) {
    for (int j = 0; j < 5; ++j) {
      outPar_reg[j] = outPar[n+j*opN]; 
    }
    for (int j = 0; j < LL; ++j) {
      /*sh_errorProp[j][i] = errorProp[n + j*epN]; */
      /*sh_errorProp[j] = errorProp[n + j*epN]; */
    }
   //const T& ConstAt(idx_t n, idx_t i, idx_t j) const { return fArray[(i * D2 + j) * N + n]; }
    const float& xin = inPar[n + 0*ipN]; 
    const float& yin = inPar[n + 1*ipN]; 
    const float& pxin = inPar[n + 3*ipN]; 
    const float& pyin = inPar[n + 4*ipN]; 
    const float& pzin = inPar[n + 5*ipN]; 
    const float& r = msRad; 
    float r0 = hipo_merged(xin, yin);

    if (fabs(r-r0)<0.0001) {
      // get an identity matrix
      computeJacobianSimple_merged(n, sh_errorProp, epN, 0, 1, 1, 1, 1, 1, 0, 1, 0, N);
      return;  // continue;
    }
    float pt2    = pxin*pxin+pyin*pyin;
    float pt     = sqrt(pt2);
    float ptinv  = 1./pt;
    float pt2inv = ptinv*ptinv;
    //p=0.3Br => r=p/(0.3*B)
    float k = inChg[n] * 100. / (-0.299792458*Config::Bfield);
    float invcurvature = 1./(pt*k);//in 1./cm
    float ctgTheta=pzin*ptinv;

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
    // TODO: remove this function, has been superseed by the register version
    // const unsigned int Niter = 5+std::round(r-r0)/2;
    for (unsigned int iter=0; iter < Config::Niter; ++iter) {
      x  = outPar_reg[0];
      y  = outPar_reg[1];
      px = outPar_reg[3];
      py = outPar_reg[4];
      r0 = hipo_merged(outPar_reg[0], outPar_reg[1]);

      totalDistance += (r-r0);
      if (Config::useTrigApprox) {  // TODO: uncomment
        sincos4_merged((r-r0)*invcurvature, sinAP, cosAP);
      } else {
        cosAP=cos((r-r0)*invcurvature);
        sinAP=sin((r-r0)*invcurvature);
      }

      //helix propagation formulas
      //http://www.phys.ufl.edu/~avery/fitting/fitting4.pdf
      outPar_reg[0] = outPar_reg[0] + k*(px*sinAP-py*(1-cosAP));
      outPar_reg[1] = outPar_reg[1] + k*(py*sinAP+px*(1-cosAP));
      outPar_reg[2] = outPar_reg[2] + (r-r0)*ctgTheta;
      outPar_reg[3] = px*cosAP-py*sinAP;
      outPar_reg[4] = py*cosAP+px*sinAP;
      //outPar.At(n, 5, 0) = pz; //take this out as it is redundant

      if (Config::useSimpleJac==0 && 
          iter +1 != Config::Niter &&
          r0 > 0 && fabs((r-r0)*invcurvature)>0.000000001) {
        //update derivatives on total distance for next step, where totalDistance+=r-r0
        //now r0 depends on px and py
        r0 = 1./r0;//WARNING, now r0 is r0inv (one less temporary)

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
      float& TD=totalDistance;
      float  TP=TD*invcurvature;//totalAngPath
      
      float& iC=invcurvature;
      float dCdpx = k*pxin*ptinv;
      float dCdpy = k*pyin*ptinv;
      float dTPdx = dTDdx*iC;
      float dTPdy = dTDdy*iC;
      float dTPdpx = (dTDdpx - TD*dCdpx*iC)*iC; // MT change: avoid division
      float dTPdpy = (dTDdpy - TD*dCdpy*iC)*iC; // MT change: avoid division
      
      float cosTP, sinTP;
      if (Config::useTrigApprox) {
        sincos4_merged(TP, sinTP, cosTP);
      } else {
        cosTP = cos(TP);
        sinTP = sin(TP);
      }

      if (Config::useSimpleJac) { 
        //assume total path length s as given and with no uncertainty
        float p = pt2 + pzin*pzin;
        p = sqrt(p);
        float s = TD*p*ptinv;
        computeJacobianSimple_merged(n, sh_errorProp, epN, s, k, p, pxin, pyin, pzin, TP, cosTP, sinTP, N);
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

        sh_errorProp[(0*L + 0)] = 1 + k*dTPdx*(pxin*cosTP - pyin*sinTP);	//dxdx;
        sh_errorProp[(0*L + 1)] = k*dTPdy*(pxin*cosTP - pyin*sinTP);	//dxdy;
        sh_errorProp[(0*L + 2)] = 0.;
        sh_errorProp[(0*L + 3)] = k*(sinTP + pxin*cosTP*dTPdpx - pyin*sinTP*dTPdpx); //dxdpx;
        sh_errorProp[(0*L + 4)] = k*(pxin*cosTP*dTPdpy - 1. + cosTP - pyin*sinTP*dTPdpy);//dxdpy;
        sh_errorProp[(0*L + 5)] = 0.;

        sh_errorProp[(1*L + 0)] = k*dTPdx*(pyin*cosTP + pxin*sinTP);	//dydx;
        sh_errorProp[(1*L + 1)] = 1 + k*dTPdy*(pyin*cosTP + pxin*sinTP);	//dydy;
        sh_errorProp[(1*L + 2)] = 0.;
        sh_errorProp[(1*L + 3)] = k*(pyin*cosTP*dTPdpx + 1. - cosTP + pxin*sinTP*dTPdpx);//dydpx;
        sh_errorProp[(1*L + 4)] = k*(sinTP + pyin*cosTP*dTPdpy + pxin*sinTP*dTPdpy); //dydpy;
        sh_errorProp[(1*L + 5)] = 0.;

        sh_errorProp[(2*L + 0)] = k*pzin*dTPdx;	//dzdx;
        sh_errorProp[(2*L + 1)] = k*pzin*dTPdy;	//dzdy;
        sh_errorProp[(2*L + 2)] = 1.;
        sh_errorProp[(2*L + 3)] = k*pzin*dTPdpx;//dzdpx;
        sh_errorProp[(2*L + 4)] = k*pzin*dTPdpy;//dzdpy;
        sh_errorProp[(2*L + 5)] = k*TP; //dzdpz;

        sh_errorProp[(3*L + 0)] = -dTPdx*(pxin*sinTP + pyin*cosTP);	//dpxdx;
        sh_errorProp[(3*L + 1)] = -dTPdy*(pxin*sinTP + pyin*cosTP);	//dpxdy;
        sh_errorProp[(3*L + 2)] = 0.;
        sh_errorProp[(3*L + 3)] = cosTP - dTPdpx*(pxin*sinTP + pyin*cosTP); //dpxdpx;
        sh_errorProp[(3*L + 4)] = -sinTP - dTPdpy*(pxin*sinTP + pyin*cosTP);//dpxdpy;
        sh_errorProp[(3*L + 5)] = 0.;

        sh_errorProp[(4*L + 0)] = -dTPdx*(pyin*sinTP - pxin*cosTP); //dpydx;
        sh_errorProp[(4*L + 1)] = -dTPdy*(pyin*sinTP - pxin*cosTP);	//dpydy;
        sh_errorProp[(4*L + 2)] = 0.;
        sh_errorProp[(4*L + 3)] = +sinTP - dTPdpx*(pyin*sinTP - pxin*cosTP);//dpydpx;
        sh_errorProp[(4*L + 4)] = +cosTP - dTPdpy*(pyin*sinTP - pxin*cosTP);//dpydpy;
        sh_errorProp[(4*L + 5)] = 0.;

        sh_errorProp[(5*L + 0)] = 0.;
        sh_errorProp[(5*L + 1)] = 0.;
        sh_errorProp[(5*L + 2)] = 0.;
        sh_errorProp[(5*L + 3)] = 0.;
        sh_errorProp[(5*L + 4)] = 0.;
        sh_errorProp[(5*L + 5)] = 1.;
      }
    }
    for (int j = 0; j < 5; ++j) {
      outPar[n + j*opN] = outPar_reg[j];
    }
    for (int j = 0; j < LL; ++j) {
      /*errorProp[n + j*epN] = sh_errorProp[j][i];*/
      /*errorProp[n + j*epN] = sh_errorProp[j];*/
    }
  }
}

/// Similarity ////////////////////////////////////////////////////////////////
__device__ void similarity_fn(float* a, size_t errorProp_stride,
    float *b, size_t stride_outErr, int N) {
  size_t aN = errorProp_stride;
  size_t bN = stride_outErr;
  
  float sh_b[LL];
  float tmp[6];

  /*int i = threadIdx.x;*/
  int n = threadIdx.x + blockIdx.x * blockDim.x;

  if (n < N) {
    for (int j = 0; j < LL; j++) {
      /*sh_a[j][i] = a[n + j*aN];*/
    }
    for (int j = 0; j < LS; j++) {
      /*sh_b[j][i] = b[n + j*bN];*/
      sh_b[j] = b[n + j*bN];
    }

    tmp[ 0] = a[0]*sh_b[ 0] + a[1]*sh_b[ 1] + a[3]*sh_b[ 6] + a[4]*sh_b[10];
    tmp[ 1] = a[0]*sh_b[ 1] + a[1]*sh_b[ 2] + a[3]*sh_b[ 7] + a[4]*sh_b[11];
    /*tmp[ 2] = a[0]*sh_b[ 3] + a[1]*sh_b[ 4] + a[3]*sh_b[ 8] + a[4]*sh_b[12];*/
    tmp[ 3] = a[0]*sh_b[ 6] + a[1]*sh_b[ 7] + a[3]*sh_b[ 9] + a[4]*sh_b[13];
    tmp[ 4] = a[0]*sh_b[10] + a[1]*sh_b[11] + a[3]*sh_b[13] + a[4]*sh_b[14];
    /*tmp[ 5] = a[0]*sh_b[15] + a[1]*sh_b[16] + a[3]*sh_b[18] + a[4]*sh_b[19];*/

    b[ 0*bN+n] = tmp[ 0]*a[0] + tmp[ 1]*a[1] + tmp[ 3]*a[3] + tmp[ 4]*a[4];


    tmp[ 0] = a[6]*sh_b[ 0] + a[7]*sh_b[ 1] + a[9]*sh_b[ 6] + a[10]*sh_b[10];
    tmp[ 1] = a[6]*sh_b[ 1] + a[7]*sh_b[ 2] + a[9]*sh_b[ 7] + a[10]*sh_b[11];
    /*tmp[ 8] = a[6]*sh_b[ 3] + a[7]*sh_b[ 4] + a[9]*sh_b[ 8] + a[10]*sh_b[12];*/
    tmp[ 3] = a[6]*sh_b[ 6] + a[7]*sh_b[ 7] + a[9]*sh_b[ 9] + a[10]*sh_b[13];
    tmp[ 4] = a[6]*sh_b[10] + a[7]*sh_b[11] + a[9]*sh_b[13] + a[10]*sh_b[14];
    /*tmp[11] = a[6]*sh_b[15] + a[7]*sh_b[16] + a[9]*sh_b[18] + a[10]*sh_b[19];*/

    b[ 1*bN+n] = tmp[ 0]*a[0] + tmp[ 1]*a[1] + tmp[ 3]*a[3] + tmp[ 4]*a[4];
    b[ 2*bN+n] = tmp[ 0]*a[6] + tmp[ 1]*a[7] + tmp[ 3]*a[9] + tmp[ 4]*a[10];


    tmp[ 0] = a[12]*sh_b[ 0] + a[13]*sh_b[ 1] + sh_b[ 3] + a[15]*sh_b[ 6] + a[16]*sh_b[10] + a[17]*sh_b[15];
    tmp[ 1] = a[12]*sh_b[ 1] + a[13]*sh_b[ 2] + sh_b[ 4] + a[15]*sh_b[ 7] + a[16]*sh_b[11] + a[17]*sh_b[16];
    tmp[ 2] = a[12]*sh_b[ 3] + a[13]*sh_b[ 4] + sh_b[ 5] + a[15]*sh_b[ 8] + a[16]*sh_b[12] + a[17]*sh_b[17];
    tmp[ 3] = a[12]*sh_b[ 6] + a[13]*sh_b[ 7] + sh_b[ 8] + a[15]*sh_b[ 9] + a[16]*sh_b[13] + a[17]*sh_b[18];
    tmp[ 4] = a[12]*sh_b[10] + a[13]*sh_b[11] + sh_b[12] + a[15]*sh_b[13] + a[16]*sh_b[14] + a[17]*sh_b[19];
    tmp[ 5] = a[12]*sh_b[15] + a[13]*sh_b[16] + sh_b[17] + a[15]*sh_b[18] + a[16]*sh_b[19] + a[17]*sh_b[20];

    b[ 3*bN+n] = tmp[ 0]*a[0] + tmp[ 1]*a[1]           + tmp[ 3]*a[3] + tmp[ 4]*a[4];
    b[ 4*bN+n] = tmp[ 0]*a[6] + tmp[ 1]*a[7]           + tmp[ 3]*a[9] + tmp[ 4]*a[10];
    b[ 5*bN+n] = tmp[ 0]*a[12] + tmp[ 1]*a[13] + tmp[ 2] + tmp[ 3]*a[15] + tmp[ 4]*a[16] + tmp[ 5]*a[17];


    tmp[ 0] = a[18]*sh_b[ 0] + a[19]*sh_b[ 1] + a[21]*sh_b[ 6] + a[22]*sh_b[10];
    tmp[ 1] = a[18]*sh_b[ 1] + a[19]*sh_b[ 2] + a[21]*sh_b[ 7] + a[22]*sh_b[11];
    tmp[ 2] = a[18]*sh_b[ 3] + a[19]*sh_b[ 4] + a[21]*sh_b[ 8] + a[22]*sh_b[12];
    tmp[ 3] = a[18]*sh_b[ 6] + a[19]*sh_b[ 7] + a[21]*sh_b[ 9] + a[22]*sh_b[13];
    tmp[ 4] = a[18]*sh_b[10] + a[19]*sh_b[11] + a[21]*sh_b[13] + a[22]*sh_b[14];
    tmp[ 5] = a[18]*sh_b[15] + a[19]*sh_b[16] + a[21]*sh_b[18] + a[22]*sh_b[19];

    b[ 6*bN+n] = tmp[ 0]*a[0] + tmp[ 1]*a[1]           + tmp[ 3]*a[3] + tmp[ 4]*a[4];
    b[ 7*bN+n] = tmp[ 0]*a[6] + tmp[ 1]*a[7]           + tmp[ 3]*a[9] + tmp[ 4]*a[10];
    b[ 8*bN+n] = tmp[ 0]*a[12] + tmp[ 1]*a[13] + tmp[ 2] + tmp[ 3]*a[15] + tmp[ 4]*a[16] + tmp[ 5]*a[17];
    b[ 9*bN+n] = tmp[ 0]*a[18] + tmp[ 1]*a[19]           + tmp[ 3]*a[21] + tmp[ 4]*a[22];


    tmp[ 0] = a[24]*sh_b[ 0] + a[25]*sh_b[ 1] + a[27]*sh_b[ 6] + a[28]*sh_b[10];
    tmp[ 1] = a[24]*sh_b[ 1] + a[25]*sh_b[ 2] + a[27]*sh_b[ 7] + a[28]*sh_b[11];
    tmp[ 2] = a[24]*sh_b[ 3] + a[25]*sh_b[ 4] + a[27]*sh_b[ 8] + a[28]*sh_b[12];
    tmp[ 3] = a[24]*sh_b[ 6] + a[25]*sh_b[ 7] + a[27]*sh_b[ 9] + a[28]*sh_b[13];
    tmp[ 4] = a[24]*sh_b[10] + a[25]*sh_b[11] + a[27]*sh_b[13] + a[28]*sh_b[14];
    tmp[ 5] = a[24]*sh_b[15] + a[25]*sh_b[16] + a[27]*sh_b[18] + a[28]*sh_b[19];

    b[10*bN+n] = tmp[ 0]*a[0] + tmp[ 1]*a[1]           + tmp[ 3]*a[3] + tmp[ 4]*a[4];
    b[11*bN+n] = tmp[ 0]*a[6] + tmp[ 1]*a[7]           + tmp[ 3]*a[9] + tmp[ 4]*a[10];
    b[12*bN+n] = tmp[ 0]*a[12] + tmp[ 1]*a[13] + tmp[ 2] + tmp[ 3]*a[15] + tmp[ 4]*a[16] + tmp[ 5]*a[17];
    b[13*bN+n] = tmp[ 0]*a[18] + tmp[ 1]*a[19]           + tmp[ 3]*a[21] + tmp[ 4]*a[22];
    b[14*bN+n] = tmp[ 0]*a[24] + tmp[ 1]*a[25]           + tmp[ 3]*a[27] + tmp[ 4]*a[28];

    tmp[ 0] = sh_b[15];
    tmp[ 1] = sh_b[16];
    tmp[ 2] = sh_b[17];
    tmp[ 3] = sh_b[18];
    tmp[ 4] = sh_b[19];
    tmp[ 5] = sh_b[20];

    // MultHelixPropTransp
    b[15*bN+n] = tmp[ 0]*a[0] + tmp[ 1]*a[1]           + tmp[ 3]*a[3] + tmp[ 4]*a[4];
    b[16*bN+n] = tmp[ 0]*a[6] + tmp[ 1]*a[7]           + tmp[ 3]*a[9] + tmp[ 4]*a[10];
    b[17*bN+n] = tmp[ 0]*a[12] + tmp[ 1]*a[13] + tmp[ 2] + tmp[ 3]*a[15] + tmp[ 4]*a[16] + tmp[ 5]*a[17];
    b[18*bN+n] = tmp[ 0]*a[18] + tmp[ 1]*a[19]           + tmp[ 3]*a[21] + tmp[ 4]*a[22];
    b[19*bN+n] = tmp[ 0]*a[24] + tmp[ 1]*a[25]           + tmp[ 3]*a[27] + tmp[ 4]*a[28];
    b[20*bN+n] = tmp[ 5];
  }
}

__global__ void propagationMerged_kernel(float *msPar, size_t stride_msPar, 
    float *inPar, size_t inPar_stride, int *inChg,
    float *outPar, size_t outPar_stride, float *errorProp,
    size_t errorProp_stride, float *outErr, size_t outErr_stride, int N) {

  int n = threadIdx.x + blockIdx.x * blockDim.x;
  float msRad_reg;
  // Using registers instead of shared memory is ~ 30% faster.
  float sh_errorProp[LL];
  if (n < N) {
    computeMsRad_fn(msPar, stride_msPar, &msRad_reg, N);
    if (Config::doIterative) {
      helixAtRFromIterative_fn(inPar, inPar_stride,
          inChg, outPar, outPar_stride, msRad_reg, 
          sh_errorProp, errorProp_stride, N);
    }
    similarity_fn(sh_errorProp, errorProp_stride, outErr, outErr_stride, N);
  }
}


void propagationMerged_wrapper(cudaStream_t& stream,
    GPlex<float>& msPar,
    GPlex<float>& inPar, GPlex<int>& inChg,
    GPlex<float>& outPar, GPlex<float>& errorProp,
    GPlex<float>& outErr, 
    const int N) {
  int gridx = std::min((N-1)/BLOCK_SIZE_X + 1,
                       MAX_BLOCKS_X);
  dim3 grid(gridx, 1, 1);
  dim3 block(BLOCK_SIZE_X, 1, 1);
  propagationMerged_kernel <<<grid, block, 0, stream >>>
    (msPar.ptr, msPar.stride,
     inPar.ptr, inPar.stride, inChg.ptr,
     outPar.ptr, outPar.stride, errorProp.ptr,
     errorProp.stride, outErr.ptr, outErr.stride, N);
}
