#include "PropagationMPlex.h"

//==============================================================================
// propagateLineToRMPlex
//==============================================================================

void propagateLineToRMPlex(const MPlexLS &psErr,  const MPlexLV& psPar,
                           const MPlexHS &msErr,  const MPlexHV& msPar,
                           MPlexLS &outErr,       MPlexLV& outPar)
{
   // XXX Regenerate parts below with a script.

   const idx_t N  = NN;

#pragma simd
   for (int n = 0; n < N; ++n)
   {
      const float dr  = hipo(msPar[0 * N + n], msPar[1 * N + n]) - hipo(psPar[0 * N + n], psPar[1 * N + n]);

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

}

// This is a tough one to crack ... sigh.
// 1. compute radius
// 2. pass in charge, too (line didn;t need it)
// 3. matriplexofy temporaries ???
// 4. all access to matriplexes by index, somehow ... sigh

// void propagateHelixToRMPlex(TrackState& inputState, float r, TrackState& result)
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

#pragma simd
   for (int n = 0; n < N; ++n)
   {
      const float& xin  = inPar.ConstAt(n, 0, 0);
      const float& yin  = inPar.ConstAt(n, 1, 0);
      const float& pxin = inPar.ConstAt(n, 3, 0);
      const float& pyin = inPar.ConstAt(n, 4, 0);
      const float& pzin = inPar.ConstAt(n, 5, 0);

      float r    = hipo(msPar.ConstAt(n, 0, 0), msPar.ConstAt(n, 1, 0));
      float r0in = hipo(xin, yin);

#ifdef DEBUG
      if (dump) std::cout << "attempt propagation from r=" << r0in << " to r=" << r << std::endl;
      if (dump) std::cout << "x=" << xin << " y=" << yin << " px=" << pxin << " py=" << pyin << " pz=" << pzin << " q=" << inChg.ConstAt(n, 0, 0) << std::endl;
      if ((r0in-r)>=0) {
         if (dump) std::cout << "target radius same or smaller than starting point, returning input" << std::endl;
         return;
      }
#endif

      float pt2    = pxin*pxin+pyin*pyin;
      float pt     = sqrt(pt2);
      float ptinv  = 1./pt;
      float pt2inv = ptinv*ptinv;
      //p=0.3Br => r=p/(0.3*B)
      float k = inChg.ConstAt(n, 0, 0) * 100. / (-0.299792458*3.8);
      float invcurvature = 1./(pt*k);//in 1./cm

#ifdef DEBUG
      if (dump) std::cout << "curvature=" << 1./invcurvature << std::endl;
#endif

      float ctgTheta=pzin*ptinv;
      //variables to be updated at each iterations
      //derivatives initialized to value for first iteration, i.e. distance = r-r0in
      float totalDistance = 0;
      float dTDdx = r0in>0. ? -xin/r0in : 0.;
      float dTDdy = r0in>0. ? -yin/r0in : 0.;
      float dTDdpx = 0.;
      float dTDdpy = 0.;
      //temporaries used within the loop (declare here to reduce memory operations)
      float x = 0.;
      float y = 0.;
      float px = 0.;
      float py = 0.;
      float r0 = 0.;
      float angPath = 0.;//maybe this temporary can be removed?
      float cosAP=0.;
      float sinAP=0.;
      float dAPdx = 0.;
      float dAPdy = 0.;
      float dAPdpx = 0.;
      float dAPdpy = 0.;
      // float dxdvar = 0.;
      // float dydvar = 0.;
      //5 iterations is a good starting point
      const unsigned int Niter = 5;
      for (unsigned int i=0;i<Niter;++i)
      {
#ifdef DEBUG
         if (dump) std::cout << "propagation iteration #" << i << std::endl;
#endif

         x  = outPar.At(n, 0, 0);
         y  = outPar.At(n, 1, 0);
         px = outPar.At(n, 3, 0);
         py = outPar.At(n, 4, 0);
         r0 = hipo(outPar.At(n, 0, 0), outPar.At(n, 1, 0));

#ifdef DEBUG
         if (dump) std::cout << "r0=" << r0 << " pt=" << pt << std::endl;
         if (dump) {
            if (r==r0) {
               std::cout << "distance = 0 at iteration=" << i << std::endl;
               break;
            }
         }
#endif

         //distance=r-r0;//remove temporary
         totalDistance+=(r-r0);
#ifdef DEBUG
         if (dump) std::cout << "distance=" << (r-r0) << std::endl;
#endif
         angPath = (r-r0)*invcurvature;
#ifdef DEBUG
         if (dump) std::cout << "angPath=" << angPath << std::endl;
#endif

         // cosAP=cos(angPath);
         // sinAP=sin(angPath);
         sincos4(angPath, sinAP, cosAP);

         //helix propagation formulas
         //http://www.phys.ufl.edu/~avery/fitting/fitting4.pdf
         outPar.At(n, 0, 0) = outPar.At(n, 0, 0) + k*(px*sinAP-py*(1-cosAP));
         outPar.At(n, 1, 0) = outPar.At(n, 1, 0) + k*(py*sinAP+px*(1-cosAP));
         outPar.At(n, 2, 0) = outPar.At(n, 2, 0) + (r-r0)*ctgTheta;
         outPar.At(n, 3, 0) = px*cosAP-py*sinAP;
         outPar.At(n, 4, 0) = py*cosAP+px*sinAP;
         //outPar.At(n, 5, 0) = pz; //take this out as it is redundant

         if (i+1 != Niter && r0 > 0)
         {
            //update derivatives on total distance for next step, where totalDistance+=r-r0
            //now r0 depends on px and py
            r0 = 1./r0;//WARNING, now r0 is r0inv (one less temporary)

#ifdef DEBUG
            if (dump) std::cout << "r0=" << 1./r0 << " r0inv=" << r0 << " pt=" << pt << std::endl;
#endif

            //update derivative on D
            dAPdx = -x*r0*invcurvature;
            dAPdy = -y*r0*invcurvature;
            dAPdpx = -angPath*px*pt2inv;
            dAPdpy = -angPath*py*pt2inv;
            //reduce temporary variables
            //dxdx = 1 + k*dAPdx*(px*sinAP + py*cosAP);
            //dydx = k*dAPdx*(py*sinAP - px*cosAP);
            //dTDdx -= r0*(x*dxdx + y*dydx);
            dTDdx -= r0*(x*(1 + k*dAPdx*(px*sinAP + py*cosAP)) + y*(k*dAPdx*(py*sinAP - px*cosAP)));
            //reuse same temporary variables
            //dxdy = k*dAPdy*(px*sinAP + py*cosAP);
            //dydy = 1 + k*dAPdy*(py*sinAP - px*cosAP);
            //dTDdy -= r0*(x*dxdy + y*dydy);
            dTDdy -= r0*(x*(k*dAPdy*(px*sinAP + py*cosAP)) + y*(1 + k*dAPdy*(py*sinAP - px*cosAP)));
            //dxdpx = k*(sinAP + px*cosAP*dAPdpx - py*sinAP*dAPdpx);
            //dydpx = k*(py*cosAP*dAPdpx + 1. - cosAP + px*sinAP*dAPdpx);
            //dTDdpx -= r0*(x*dxdpx + y*dTDdpx);
            dTDdpx -= r0*(x*(k*(sinAP + px*cosAP*dAPdpx - py*sinAP*dAPdpx)) + y*(k*(py*cosAP*dAPdpx + 1. - cosAP + px*sinAP*dAPdpx)));
            //dxdpy = k*(px*cosAP*dAPdpy - 1. + cosAP - py*sinAP*dAPdpy);
            //dydpy = k*(sinAP + py*cosAP*dAPdpy + px*sinAP*dAPdpy);
            //dTDdpy -= r0*(x*dxdpy + y*(k*dydpy);
            dTDdpy -= r0*(x*(k*(px*cosAP*dAPdpy - 1. + cosAP - py*sinAP*dAPdpy)) + y*(k*(sinAP + py*cosAP*dAPdpy + px*sinAP*dAPdpy)));
         }

#ifdef DEBUG
         if (dump) std::cout << outPar.At(n, 0, 0) << " " << outPar.At(n, 1, 0) << " " << outPar.At(n, 2, 0) << std::endl;
         if (dump) std::cout << outPar.At(n, 3, 0) << " " << outPar.At(n, 4, 0) << " " << outPar.At(n, 5, 0) << std::endl;
#endif
      }

      float totalAngPath=totalDistance*invcurvature;
      float& TD=totalDistance;
      float& TP=totalAngPath;
      float& iC=invcurvature;

#ifdef DEBUG
      if (dump) std::cout << "TD=" << TD << " TP=" << TP << " arrived at r=" << sqrt(outPar.At(n, 0, 0)*outPar.At(n, 0, 0)+outPar.At(n, 1, 0)*outPar.At(n, 1, 0)) << std::endl;
#endif

      float dCdpx = k*pxin*ptinv;
      float dCdpy = k*pyin*ptinv;
      float dTPdx = dTDdx*iC;
      float dTPdy = dTDdy*iC;
      float dTPdpx = (dTDdpx - TD*dCdpx*iC)*iC; // MT change: avoid division
      float dTPdpy = (dTDdpy - TD*dCdpy*iC)*iC; // MT change: avoid division

      // float cosTP = cos(TP);
      // float sinTP = sin(TP);
      float cosTP, sinTP;
      sincos4(TP, sinTP, cosTP);

      //derive these to compute jacobian
      //x = xin + k*(pxin*sinTP-pyin*(1-cosTP));
      //y = yin + k*(pyin*sinTP+pxin*(1-cosTP));
      //z = zin + TD*ctgTheta;
      //px = pxin*cosTP-pyin*sinTP;
      //py = pyin*cosTP+pxin*sinTP;
      //pz = pzin;
      //jacobian


      // Not everything is set! GenMul patterns are used to set 0 and 1 elements.

      errorProp(n,0,0) = 1 + k*dTPdx*(pxin*sinTP + pyin*cosTP);	//dxdx;
      errorProp(n,0,1) = k*dTPdy*(pxin*sinTP + pyin*cosTP);	//dxdy;
      errorProp(n,0,3) = k*(sinTP + pxin*cosTP*dTPdpx - pyin*sinTP*dTPdpx); //dxdpx;
      errorProp(n,0,4) = k*(pxin*cosTP*dTPdpy - 1. + cosTP - pyin*sinTP*dTPdpy);//dxdpy;
      errorProp(n,1,0) = k*dTPdx*(pyin*sinTP - pxin*cosTP);	//dydx;
      errorProp(n,1,1) = 1 + k*dTPdy*(pyin*sinTP - pxin*cosTP);	//dydy;
      errorProp(n,1,3) = k*(pyin*cosTP*dTPdpx + 1. - cosTP + pxin*sinTP*dTPdpx);//dydpx;
      errorProp(n,1,4) = k*(sinTP + pyin*cosTP*dTPdpy + pxin*sinTP*dTPdpy); //dydpy;
      errorProp(n,2,0) = dTDdx*ctgTheta;	//dzdx;
      errorProp(n,2,1) = dTDdy*ctgTheta;	//dzdy;
      errorProp(n,2,3) = dTDdpx*ctgTheta - TD*pzin*pxin*pt2inv*ptinv;//dzdpx;
      errorProp(n,2,4) = dTDdpy*ctgTheta - TD*pzin*pyin*pt2inv*ptinv;//dzdpy;
      errorProp(n,2,5) = TD*ptinv; //dzdpz;
      errorProp(n,3,0) = -dTPdx*(pxin*sinTP + pyin*cosTP);	//dpxdx;
      errorProp(n,3,1) = -dTPdy*(pxin*sinTP + pyin*cosTP);	//dpxdy;
      errorProp(n,3,3) = cosTP - dTPdpx*(pxin*sinTP + pyin*cosTP); //dpxdpx;
      errorProp(n,3,4) = -sinTP - dTPdpy*(pxin*sinTP + pyin*cosTP);//dpxdpy;
      errorProp(n,4,0) = -dTPdx*(pyin*sinTP - pxin*cosTP); //dpydx;
      errorProp(n,4,1) = -dTPdy*(pyin*sinTP - pxin*cosTP);	//dpydy;
      errorProp(n,4,3) = +sinTP - dTPdpx*(pyin*sinTP - pxin*cosTP);//dpydpx;
      errorProp(n,4,4) = +cosTP - dTPdpy*(pyin*sinTP - pxin*cosTP);//dpydpy;
   }

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

   /*
     if (fabs(sqrt(outPar[0]*outPar[0]+outPar[1]*outPar[1])-r)>0.0001) {
     std::cout << "DID NOT GET TO R, dR=" << fabs(sqrt(outPar[0]*outPar[0]+outPar[1]*outPar[1])-r)
     << " r=" << r << " r0in=" << r0in << " rout=" << sqrt(outPar[0]*outPar[0]+outPar[1]*outPar[1]) << std::endl;
     std::cout << "pt=" << pt << " pz=" << inPar.At(n, 2) << std::endl;
     }
   */
}

void propagateHelixToRMPlex(const MPlexLS &inErr,  const MPlexLV& inPar,
                            const MPlexQI &inChg,  const float r,
                                  MPlexLS &outErr,       MPlexLV& outPar)
{
#ifdef DEBUG
   const bool dump = false;
#endif

   const idx_t N  = NN;

   outErr = inErr;
   outPar = inPar;

   MPlexLL errorProp;

#pragma simd
   for (int n = 0; n < N; ++n)
   {
      const float& xin  = inPar.ConstAt(n, 0, 0);
      const float& yin  = inPar.ConstAt(n, 1, 0);
      const float& pxin = inPar.ConstAt(n, 3, 0);
      const float& pyin = inPar.ConstAt(n, 4, 0);
      const float& pzin = inPar.ConstAt(n, 5, 0);

      std::cout << "propagate track from x=" << xin << " y=" << yin << " using vecunit=" << n << std::endl;

      float r0in = hipo(xin, yin);

#ifdef DEBUG
      if (dump) std::cout << "attempt propagation from r=" << r0in << " to r=" << r << std::endl;
      if (dump) std::cout << "x=" << xin << " y=" << yin << " px=" << pxin << " py=" << pyin << " pz=" << pzin << " q=" << inChg.ConstAt(n, 0, 0) << std::endl;
      if ((r0in-r)>=0) {
         if (dump) std::cout << "target radius same or smaller than starting point, returning input" << std::endl;
         return;
      }
#endif

      float pt2    = pxin*pxin+pyin*pyin;
      float pt     = sqrt(pt2);
      float ptinv  = 1./pt;
      float pt2inv = ptinv*ptinv;
      //p=0.3Br => r=p/(0.3*B)
      float k = inChg.ConstAt(n, 0, 0) * 100. / (-0.299792458*3.8);
      float invcurvature = 1./(pt*k);//in 1./cm

#ifdef DEBUG
      if (dump) std::cout << "curvature=" << 1./invcurvature << std::endl;
#endif

      float ctgTheta=pzin*ptinv;
      //variables to be updated at each iterations
      //derivatives initialized to value for first iteration, i.e. distance = r-r0in
      float totalDistance = 0;
      float dTDdx = r0in>0. ? -xin/r0in : 0.;
      float dTDdy = r0in>0. ? -yin/r0in : 0.;
      float dTDdpx = 0.;
      float dTDdpy = 0.;
      //temporaries used within the loop (declare here to reduce memory operations)
      float x = 0.;
      float y = 0.;
      float px = 0.;
      float py = 0.;
      float r0 = 0.;
      float angPath = 0.;//maybe this temporary can be removed?
      float cosAP=0.;
      float sinAP=0.;
      float dAPdx = 0.;
      float dAPdy = 0.;
      float dAPdpx = 0.;
      float dAPdpy = 0.;
      // float dxdvar = 0.;
      // float dydvar = 0.;
      //5 iterations is a good starting point
      const unsigned int Niter = 5;
      for (unsigned int i=0;i<Niter;++i)
      {
#ifdef DEBUG
         if (dump) std::cout << "propagation iteration #" << i << std::endl;
#endif

         x  = outPar.At(n, 0, 0);
         y  = outPar.At(n, 1, 0);
         px = outPar.At(n, 3, 0);
         py = outPar.At(n, 4, 0);
         r0 = hipo(outPar.At(n, 0, 0), outPar.At(n, 1, 0));

#ifdef DEBUG
         if (dump) std::cout << "r0=" << r0 << " pt=" << pt << std::endl;
         if (dump) {
            if (r==r0) {
               std::cout << "distance = 0 at iteration=" << i << std::endl;
               break;
            }
         }
#endif

         //distance=r-r0;//remove temporary
         totalDistance+=(r-r0);
#ifdef DEBUG
         if (dump) std::cout << "distance=" << (r-r0) << std::endl;
#endif
         angPath = (r-r0)*invcurvature;
#ifdef DEBUG
         if (dump) std::cout << "angPath=" << angPath << std::endl;
#endif

         // cosAP=cos(angPath);
         // sinAP=sin(angPath);
         sincos4(angPath, sinAP, cosAP);

         //helix propagation formulas
         //http://www.phys.ufl.edu/~avery/fitting/fitting4.pdf
         outPar.At(n, 0, 0) = outPar.At(n, 0, 0) + k*(px*sinAP-py*(1-cosAP));
         outPar.At(n, 1, 0) = outPar.At(n, 1, 0) + k*(py*sinAP+px*(1-cosAP));
         outPar.At(n, 2, 0) = outPar.At(n, 2, 0) + (r-r0)*ctgTheta;
         outPar.At(n, 3, 0) = px*cosAP-py*sinAP;
         outPar.At(n, 4, 0) = py*cosAP+px*sinAP;
         //outPar.At(n, 5, 0) = pz; //take this out as it is redundant

         if (i+1 != Niter && r0 > 0)
         {
            //update derivatives on total distance for next step, where totalDistance+=r-r0
            //now r0 depends on px and py
            r0 = 1./r0;//WARNING, now r0 is r0inv (one less temporary)

#ifdef DEBUG
            if (dump) std::cout << "r0=" << 1./r0 << " r0inv=" << r0 << " pt=" << pt << std::endl;
#endif

            //update derivative on D
            dAPdx = -x*r0*invcurvature;
            dAPdy = -y*r0*invcurvature;
            dAPdpx = -angPath*px*pt2inv;
            dAPdpy = -angPath*py*pt2inv;
            //reduce temporary variables
            //dxdx = 1 + k*dAPdx*(px*sinAP + py*cosAP);
            //dydx = k*dAPdx*(py*sinAP - px*cosAP);
            //dTDdx -= r0*(x*dxdx + y*dydx);
            dTDdx -= r0*(x*(1 + k*dAPdx*(px*sinAP + py*cosAP)) + y*(k*dAPdx*(py*sinAP - px*cosAP)));
            //reuse same temporary variables
            //dxdy = k*dAPdy*(px*sinAP + py*cosAP);
            //dydy = 1 + k*dAPdy*(py*sinAP - px*cosAP);
            //dTDdy -= r0*(x*dxdy + y*dydy);
            dTDdy -= r0*(x*(k*dAPdy*(px*sinAP + py*cosAP)) + y*(1 + k*dAPdy*(py*sinAP - px*cosAP)));
            //dxdpx = k*(sinAP + px*cosAP*dAPdpx - py*sinAP*dAPdpx);
            //dydpx = k*(py*cosAP*dAPdpx + 1. - cosAP + px*sinAP*dAPdpx);
            //dTDdpx -= r0*(x*dxdpx + y*dTDdpx);
            dTDdpx -= r0*(x*(k*(sinAP + px*cosAP*dAPdpx - py*sinAP*dAPdpx)) + y*(k*(py*cosAP*dAPdpx + 1. - cosAP + px*sinAP*dAPdpx)));
            //dxdpy = k*(px*cosAP*dAPdpy - 1. + cosAP - py*sinAP*dAPdpy);
            //dydpy = k*(sinAP + py*cosAP*dAPdpy + px*sinAP*dAPdpy);
            //dTDdpy -= r0*(x*dxdpy + y*(k*dydpy);
            dTDdpy -= r0*(x*(k*(px*cosAP*dAPdpy - 1. + cosAP - py*sinAP*dAPdpy)) + y*(k*(sinAP + py*cosAP*dAPdpy + px*sinAP*dAPdpy)));
         }

#ifdef DEBUG
         if (dump) std::cout << outPar.At(n, 0, 0) << " " << outPar.At(n, 1, 0) << " " << outPar.At(n, 2, 0) << std::endl;
         if (dump) std::cout << outPar.At(n, 3, 0) << " " << outPar.At(n, 4, 0) << " " << outPar.At(n, 5, 0) << std::endl;
#endif
      }

      float totalAngPath=totalDistance*invcurvature;
      float& TD=totalDistance;
      float& TP=totalAngPath;
      float& iC=invcurvature;

#ifdef DEBUG
      if (dump) std::cout << "TD=" << TD << " TP=" << TP << " arrived at r=" << sqrt(outPar.At(n, 0, 0)*outPar.At(n, 0, 0)+outPar.At(n, 1, 0)*outPar.At(n, 1, 0)) << std::endl;
#endif

      float dCdpx = k*pxin*ptinv;
      float dCdpy = k*pyin*ptinv;
      float dTPdx = dTDdx*iC;
      float dTPdy = dTDdy*iC;
      float dTPdpx = (dTDdpx - TD*dCdpx*iC)*iC; // MT change: avoid division
      float dTPdpy = (dTDdpy - TD*dCdpy*iC)*iC; // MT change: avoid division

      // float cosTP = cos(TP);
      // float sinTP = sin(TP);
      float cosTP, sinTP;
      sincos4(TP, sinTP, cosTP);

      //derive these to compute jacobian
      //x = xin + k*(pxin*sinTP-pyin*(1-cosTP));
      //y = yin + k*(pyin*sinTP+pxin*(1-cosTP));
      //z = zin + TD*ctgTheta;
      //px = pxin*cosTP-pyin*sinTP;
      //py = pyin*cosTP+pxin*sinTP;
      //pz = pzin;
      //jacobian


      // Not everything is set! GenMul patterns are used to set 0 and 1 elements.

      errorProp(n,0,0) = 1 + k*dTPdx*(pxin*sinTP + pyin*cosTP);	//dxdx;
      errorProp(n,0,1) = k*dTPdy*(pxin*sinTP + pyin*cosTP);	//dxdy;
      errorProp(n,0,3) = k*(sinTP + pxin*cosTP*dTPdpx - pyin*sinTP*dTPdpx); //dxdpx;
      errorProp(n,0,4) = k*(pxin*cosTP*dTPdpy - 1. + cosTP - pyin*sinTP*dTPdpy);//dxdpy;
      errorProp(n,1,0) = k*dTPdx*(pyin*sinTP - pxin*cosTP);	//dydx;
      errorProp(n,1,1) = 1 + k*dTPdy*(pyin*sinTP - pxin*cosTP);	//dydy;
      errorProp(n,1,3) = k*(pyin*cosTP*dTPdpx + 1. - cosTP + pxin*sinTP*dTPdpx);//dydpx;
      errorProp(n,1,4) = k*(sinTP + pyin*cosTP*dTPdpy + pxin*sinTP*dTPdpy); //dydpy;
      errorProp(n,2,0) = dTDdx*ctgTheta;	//dzdx;
      errorProp(n,2,1) = dTDdy*ctgTheta;	//dzdy;
      errorProp(n,2,3) = dTDdpx*ctgTheta - TD*pzin*pxin*pt2inv*ptinv;//dzdpx;
      errorProp(n,2,4) = dTDdpy*ctgTheta - TD*pzin*pyin*pt2inv*ptinv;//dzdpy;
      errorProp(n,2,5) = TD*ptinv; //dzdpz;
      errorProp(n,3,0) = -dTPdx*(pxin*sinTP + pyin*cosTP);	//dpxdx;
      errorProp(n,3,1) = -dTPdy*(pxin*sinTP + pyin*cosTP);	//dpxdy;
      errorProp(n,3,3) = cosTP - dTPdpx*(pxin*sinTP + pyin*cosTP); //dpxdpx;
      errorProp(n,3,4) = -sinTP - dTPdpy*(pxin*sinTP + pyin*cosTP);//dpxdpy;
      errorProp(n,4,0) = -dTPdx*(pyin*sinTP - pxin*cosTP); //dpydx;
      errorProp(n,4,1) = -dTPdy*(pyin*sinTP - pxin*cosTP);	//dpydy;
      errorProp(n,4,3) = +sinTP - dTPdpx*(pyin*sinTP - pxin*cosTP);//dpydpx;
      errorProp(n,4,4) = +cosTP - dTPdpy*(pyin*sinTP - pxin*cosTP);//dpydpy;
   }

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

   /*
     if (fabs(sqrt(outPar[0]*outPar[0]+outPar[1]*outPar[1])-r)>0.0001) {
     std::cout << "DID NOT GET TO R, dR=" << fabs(sqrt(outPar[0]*outPar[0]+outPar[1]*outPar[1])-r)
     << " r=" << r << " r0in=" << r0in << " rout=" << sqrt(outPar[0]*outPar[0]+outPar[1]*outPar[1]) << std::endl;
     std::cout << "pt=" << pt << " pz=" << inPar.At(n, 2) << std::endl;
     }
   */
}
