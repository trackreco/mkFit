#ifndef _kalmanutils_
#define _kalmanutils_

#include "Track.h"
#include "Matrix.h"

TrackState updateParameters(TrackState& propagatedState, MeasurementState& measurementState, 
			    SMatrix36& projMatrix,SMatrix63& projMatrixT) {

  bool print = false;

  SMatrix66 propErr = propagatedState.errors;
  if (print) {
    std::cout << "propErr" << std::endl;
    dumpMatrix(propErr);
  }
  SMatrix33 propErr33 = projMatrix*propagatedState.errors*projMatrixT;
  if (print) {
    std::cout << "propErr33" << std::endl;
    dumpMatrix(propErr33);
  }  

  SVector3 residual = measurementState.parameters-projMatrix*propagatedState.parameters;
  if (print) {
  std::cout << "residual: " << residual[0] << " " << residual[1] << " " << residual[2] << " "
	    << residual[3] << " " << residual[4] << " " << residual[5] << " " << std::endl;
  }
  SMatrix33 resErr = measurementState.errors+propErr33;
  if (print) {
    std::cout << "resErr" << std::endl;
    dumpMatrix(resErr);
  }

  SMatrix33 resErrInv = resErr;
  resErrInv.Invert();
  if (print) {
    std::cout << "resErrInv" << std::endl;
    dumpMatrix(resErrInv);
  }

  SMatrix63 kalmanGain = propagatedState.errors*projMatrixT*(resErrInv);
  if (print) {
    std::cout << "kalmanGain" << std::endl;
    dumpMatrix(kalmanGain);
  }
  SVector6 updatedParams = propagatedState.parameters + kalmanGain*residual;
  if (print) {
  std::cout << "updatedParams: " << updatedParams[0] << " " << updatedParams[1] << " " << updatedParams[2] << " "
	    << updatedParams[3] << " " << updatedParams[4] << " " << updatedParams[5] << " " << std::endl;
  }
  SMatrixSym66 identity = ROOT::Math::SMatrixIdentity();
  SMatrix66 updatedErrs = (identity - kalmanGain*projMatrix)*propagatedState.errors;
  if (print) {
    std::cout << "updatedErrs" << std::endl;
    dumpMatrix(updatedErrs);
  }

  TrackState result;
  result.parameters=updatedParams;
  result.errors=updatedErrs;
  return result;
}

//temporary...
void setupTrackByHand(Point& pos, Vector& mom, SMatrixSym66& covtrk, std::vector<Hit>& hits, int vers) {
  if (vers==0) {
    pos=Point(0.556433,-4.13765,1.17947);
    mom=Vector(6.93565,-99.8987,106.122);
    covtrk = ROOT::Math::SMatrixIdentity();
    covtrk(0,0)=0.00668929   ;
    covtrk(0,1)=-0.000405335 ;
    covtrk(0,2)=-0.00081875  ;
    covtrk(0,3)=-0.654305    ;
    covtrk(0,4)=8.77271      ;
    covtrk(0,5)=-9.28359     ;
    covtrk(1,1)=0.00585989   ;
    covtrk(1,2)=0.00554276   ;
    covtrk(1,3)=-0.0221493   ;
    covtrk(1,4)=0.336286     ;
    covtrk(1,5)=-0.533675    ;
    covtrk(2,2)=0.00527125   ;
    covtrk(2,3)=0.0219121    ;
    covtrk(2,4)=-0.25678     ;
    covtrk(2,5)=0.104355     ;
    covtrk(3,3)=365.307      ;
    covtrk(3,4)=-5081.05     ;
    covtrk(3,5)=5397.81      ;
    covtrk(4,4)=70741.4      ;
    covtrk(4,5)=-75154.2     ;
    covtrk(5,5)=79853.4      ; 
    Point x1(0.555814,-4.13765,1.18106);
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=7.89828e-07; 
    covx1(0,1)=1.71343e-10;
    covx1(0,2)=-3.65348e-11;
    covx1(1,1)=3.75336e-14;
    covx1(1,2)=3.71905e-11;
    covx1(2,2)=3.81485e-06;
    Hit hit1(x1,covx1);
    Point x2(0.759732,-7.0482,4.26977);
    SMatrixSym33 covx2 = ROOT::Math::SMatrixIdentity();
    covx2(0,0)=1.16538e-06; 
    covx2(0,1)=3.15313e-10;
    covx2(0,2)=-9.24943e-12;
    covx2(1,1)=8.59631e-14;
    covx2(1,2)=-4.24538e-11;
    covx2(2,2)=2.7723e-06;
    Hit hit2(x2,covx2);
    hits.push_back(hit1);
    hits.push_back(hit2);
  } else {
    pos=Point(-3.97545,-1.39023,-3.20927);
    mom=Vector(-822.213,-348.139,-765.718);
    covtrk = ROOT::Math::SMatrixIdentity();
    covtrk(0,0)=0.00339944  ;
    covtrk(0,1)=-0.000129504;
    covtrk(0,2)=-0.00359138 ;
    covtrk(0,3)=-86.3911    ;
    covtrk(0,4)=-36.5052    ;
    covtrk(0,5)=-79.7691    ;
    covtrk(1,1)=0.00386787 ;
    covtrk(1,2)=-0.0016195 ;
    covtrk(1,3)=154.307    ;
    covtrk(1,4)=65.1894    ;
    covtrk(1,5)=143.913    ;
    covtrk(2,2)=0.00459267;
    covtrk(2,2)=22.6083   ;
    covtrk(2,2)=9.55969   ;
    covtrk(2,2)=20.2235   ;
    covtrk(3,3)=4.59221e+07;
    covtrk(3,3)=1.94154e+07;
    covtrk(3,3)=4.27286e+07;
    covtrk(4,4)=8.20865e+06;
    covtrk(4,4)=1.80653e+07;
    covtrk(5,5)=3.9758e+07;
    Point x1(-3.97592,-1.38942,-3.21347);
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=2.8478e-07  ;
    covx1(0,1)=-4.93169e-07;
    covx1(0,2)=-1.68421e-11;
    covx1(1,1)=8.54048e-07 ;
    covx1(1,2)=-2.08965e-11;
    covx1(2,2)=2.82947e-06 ;
    Hit hit1(x1,covx1);
    Point x2(-4.40679,-1.57342,-3.60635);
    SMatrixSym33 covx2 = ROOT::Math::SMatrixIdentity();
    covx2(0,0)=4.19391e-09 ;
    covx2(0,1)=-2.38012e-08;
    covx2(0,2)=3.16473e-12 ;
    covx2(1,1)=1.35076e-07 ;
    covx2(1,2)=-2.23185e-11;
    covx2(2,2)=3.05883e-06 ;
    Hit hit2(x2,covx2);
    Point x3(-7.03704,-2.68632,-6.06063);
    SMatrixSym33 covx3 = ROOT::Math::SMatrixIdentity();
    covx3(0,0)=1.26699e-07 ;
    covx3(0,1)=-3.89985e-07;
    covx3(0,2)=2.77271e-12 ;
    covx3(1,1)=1.20039e-06 ;
    covx3(1,2)=-5.58307e-12;
    covx3(2,2)=2.93281e-06 ;
    Hit hit3(x3,covx3);
    Point x4(-9.67263,-3.80139,-8.51372);
    SMatrixSym33 covx4 = ROOT::Math::SMatrixIdentity();
    covx4(0,0)=2.05148e-06 ;
    covx4(0,1)=-5.22726e-06;
    covx4(0,2)=-1.94228e-10;
    covx4(1,1)=1.33193e-05 ;
    covx4(1,2)=-1.02983e-10;
    covx4(2,2)=1.88032e-05 ;
    Hit hit4(x4,covx4);
    Point x5(-21.9214,-8.98664,-24.2565);
    SMatrixSym33 covx5 = ROOT::Math::SMatrixIdentity();
    covx5(0,0)=2.1603e+28 ;
    covx5(0,1)=4.64935e+26;
    covx5(0,2)=2.71129e+33;
    covx5(1,1)=1.00062e+25;
    covx5(1,2)=5.83519e+31;
    covx5(2,2)=3.40282e+38;
    Hit hit5(x5,covx5);
    Point x6(-22.3404,-8.71725,-24.2776);
    SMatrixSym33 covx6 = ROOT::Math::SMatrixIdentity();
    covx6(0,0)=8.35807e+35;
    covx6(0,1)=-1.4616e+36;
    covx6(0,2)=1.67802e+37;
    covx6(1,1)=2.55593e+36;
    covx6(1,2)=-2.9344e+37;
    covx6(2,2)=3.36891e+38;
    Hit hit6(x6,covx6);
    hits.push_back(hit1);
    hits.push_back(hit2);
    hits.push_back(hit3);
    hits.push_back(hit4);
    hits.push_back(hit5);
    hits.push_back(hit6);
  }  
}

#endif
