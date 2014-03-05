#ifndef _propagation_
#define _propagation_

#include "Track.h"
#include "Matrix.h"

TrackState propagateLineToR(TrackState& inputState, float r) {

  bool dump = false;

  SVector6& par = inputState.parameters;
  SMatrix66& err = inputState.errors;

  //straight line for now
  float r0 = sqrt(par.At(0)*par.At(0)+par.At(1)*par.At(1));
  float dr = r-r0;
  float pt = sqrt(par.At(3)*par.At(3)+par.At(4)*par.At(4));
  float path = dr/pt;//this works only if direction is along radius, i.e. origin is at (0,0)

  TrackState result;

  SMatrix66 propMatrix = ROOT::Math::SMatrixIdentity();
  propMatrix(0,3)=path;
  propMatrix(1,4)=path;
  propMatrix(2,5)=path;
  result.parameters=propMatrix*par;
  if (dump) {
    //test R of propagation
    std::cout << "initial R=" << r0 << std::endl;
    std::cout << "target R=" << r << std::endl;
    std::cout << "arrived at R=" << sqrt(result.parameters[0]*result.parameters[0]+result.parameters[1]*result.parameters[1]) << std::endl;
  }

  result.errors=propMatrix*err*ROOT::Math::Transpose(propMatrix);  
  return result;
}

//fixme what about charge????
TrackState propagateHelixToR(TrackState& inputState, int& charge, float r) {

  bool dump = false;

  float xin = inputState.parameters.At(0);
  float yin = inputState.parameters.At(1);
  float pxin = inputState.parameters.At(3);
  float pyin = inputState.parameters.At(4);
  float pzin = inputState.parameters.At(5);

  float pt2 = pxin*pxin+pyin*pyin;
  float pt = sqrt(pt2);
  float pt3 = pt*pt2;
  //p=0.3Br => r=p/(0.3*B)
  float k=charge*100./(-0.299792458*3.8);
  float curvature = pt*k;//in cm
  if (dump) std::cout << "curvature=" << curvature << std::endl;
  float ctgTheta=pzin/pt;

  float totalDistance = 0;
  float r0in = sqrt(xin*xin+yin*yin);
  float dTDdx = r0in>0. ? -xin/r0in : 0;
  float dTDdy = r0in>0. ? -yin/r0in : 0;
  float dTDdpx = 0.;
  float dTDdpy = 0.;

  //make a copy for now...
  SVector6 par = inputState.parameters;
  SMatrix66 err = inputState.errors;

  //5 iterations is a good starting point
  unsigned int Niter = 5;
  for (unsigned int i=0;i<Niter;++i) {

    if (dump) std::cout << "propagation iteration #" << i << std::endl;

    float x = par.At(0);
    float y = par.At(1);
    float z = par.At(2);
    float px = par.At(3);
    float py = par.At(4);
    float pz = par.At(5);

    float r0 = sqrt(x*x+y*y);
    if (dump) std::cout << "r0=" << r0 << " pt=" << pt << std::endl;

    if (dump) {
      if (r==r0) {
	std::cout << "distance = 0 at iteration=" <<  i << std::endl;
	break;
      } 
    }

    float distance = r-r0;//fixme compute real distance between two points
    totalDistance+=distance;
    if (dump) std::cout << "distance=" << distance << std::endl;  
    float angPath = distance/curvature;
    if (dump) std::cout << "angPath=" << angPath << std::endl;
    float cosAP=cos(angPath);
    float sinAP=sin(angPath);

    //http://www.phys.ufl.edu/~avery/fitting/fitting4.pdf
    par.At(0) = x + k*(px*sinAP-py*(1-cosAP));
    par.At(1) = y + k*(py*sinAP+px*(1-cosAP));
    par.At(2) = z + distance*ctgTheta;

    par.At(3) = px*cosAP-py*sinAP;
    par.At(4) = py*cosAP+px*sinAP;
    par.At(5) = pz;
    
    if (i+1!=Niter && r0>0.) {
      float r0inv = 1./r0;
      if (dump) std::cout << "r0=" << r0 << " r0inv=" << r0inv << " pt=" << pt << std::endl;
      //update derivative on D
      float dAPdx = -x/(r0*curvature);
      float dAPdy = -y/(r0*curvature);
      float dAPdpx = -angPath*px/pt2;
      float dAPdpy = -angPath*py/pt2;

      float dxdx = 1 + k*dAPdx*(px*sinAP + py*cosAP);
      float dxdy = k*dAPdy*(px*sinAP + py*cosAP);
      float dydx = k*dAPdx*(py*sinAP - px*cosAP);
      float dydy = 1 + k*dAPdy*(py*sinAP - px*cosAP);

      float dxdpx = k*(sinAP + px*cosAP*dAPdpx - py*sinAP*dAPdpx);
      float dxdpy = k*(px*cosAP*dAPdpy - 1. + cosAP - py*sinAP*dAPdpy);
      float dydpx = k*(py*cosAP*dAPdpx + 1. - cosAP + px*sinAP*dAPdpx);
      float dydpy = k*(sinAP + py*cosAP*dAPdpy + px*sinAP*dAPdpy);

      dTDdx -= r0inv*(x*dxdx + y*dydx);
      dTDdy -= r0inv*(x*dxdy + y*dydy);
      dTDdpx -= r0inv*(x*dxdpx + y*dydpx);
      dTDdpy -= r0inv*(x*dxdpy + y*dydpy);
    }

    if (dump) std::cout << par.At(0) << " " << par.At(1) << " " << par.At(2) << std::endl;
    if (dump) std::cout << par.At(3) << " " << par.At(4) << " " << par.At(5) << std::endl;

  }

  float totalAngPath=totalDistance/curvature;

  float TD=totalDistance;
  float TP=totalAngPath;
  float C=curvature;

  if (dump) std::cout << "TD=" << TD << " TP=" << TP << " arrived at r=" << sqrt(par.At(0)*par.At(0)+par.At(1)*par.At(1)) << std::endl;

  float dCdpx = k*pxin/pt;
  float dCdpy = k*pyin/pt;

  float dTPdx = dTDdx/C;
  float dTPdy = dTDdy/C;
  float dTPdpx = (dTDdpx*C - TD*dCdpx)/(C*C);
  float dTPdpy = (dTDdpy*C - TD*dCdpy)/(C*C);

  float cosTP = cos(TP);
  float sinTP = sin(TP);

  //par.At(0) = xin + k*(pxin*sinTP-pyin*(1-cosTP));
  //par.At(1) = yin + k*(pyin*sinTP+pxin*(1-cosTP));
  //par.At(2) = zin + TD*ctgTheta;

  float dxdx = 1 + k*dTPdx*(pxin*sinTP + pyin*cosTP);
  float dxdy = k*dTPdy*(pxin*sinTP + pyin*cosTP);
  float dydx = k*dTPdx*(pyin*sinTP - pxin*cosTP);
  float dydy = 1 + k*dTPdy*(pyin*sinTP - pxin*cosTP);

  float dxdpx = k*(sinTP + pxin*cosTP*dTPdpx - pyin*sinTP*dTPdpx);
  float dxdpy = k*(pxin*cosTP*dTPdpy - 1. + cosTP - pyin*sinTP*dTPdpy);
  float dydpx = k*(pyin*cosTP*dTPdpx + 1. - cosTP + pxin*sinTP*dTPdpx);
  float dydpy = k*(sinTP + pyin*cosTP*dTPdpy + pxin*sinTP*dTPdpy);

  float dzdx = dTDdx*ctgTheta;
  float dzdy = dTDdy*ctgTheta;

  float dzdpx = dTDdpx*ctgTheta - TD*pzin*pxin/pt3;
  float dzdpy = dTDdpy*ctgTheta - TD*pzin*pyin/pt3;
  float dzdpz = TD/pt;//fixme if I set this term to 0 then it works...

  //par.At(3) = pxin*cosTP-pyin*sinTP;
  //par.At(4) = pyin*cosTP+pxin*sinTP;
  //par.At(5) = pzin;

  float dpxdx = -dTPdx*(pxin*sinTP + pyin*cosTP);
  float dpxdy = -dTPdy*(pxin*sinTP + pyin*cosTP);
  float dpydx = -dTPdx*(pyin*sinTP - pxin*cosTP);
  float dpydy = -dTPdy*(pyin*sinTP - pxin*cosTP);
  
  float dpxdpx = cosTP - dTPdpx*(pxin*sinTP + pyin*cosTP);
  float dpxdpy = -sinTP - dTPdpy*(pxin*sinTP + pyin*cosTP);
  float dpydpx = +sinTP - dTPdpx*(pyin*sinTP - pxin*cosTP);
  float dpydpy = cosTP - dTPdpy*(pyin*sinTP - pxin*cosTP);

  //jacobian
  SMatrix66 errorProp = ROOT::Math::SMatrixIdentity();
  errorProp(0,0)=dxdx;
  errorProp(0,1)=dxdy;
  errorProp(0,3)=dxdpx;
  errorProp(0,4)=dxdpy;

  errorProp(1,0)=dydx;
  errorProp(1,1)=dydy;
  errorProp(1,3)=dydpx;
  errorProp(1,4)=dydpy;

  errorProp(2,0)=dzdx;
  errorProp(2,1)=dzdy;
  errorProp(2,3)=dzdpx;
  errorProp(2,4)=dzdpy;
  errorProp(2,5)=dzdpz;

  errorProp(3,0)=dpxdx;
  errorProp(3,1)=dpxdy;
  errorProp(3,3)=dpxdpx;
  errorProp(3,4)=dpxdpy;

  errorProp(4,0)=dpydx;
  errorProp(4,1)=dpydy;
  errorProp(4,3)=dpydpx;
  errorProp(4,4)=dpydpy;

  if (dump) {
    std::cout << "errorProp" << std::endl;
    dumpMatrix(errorProp);
  }

  TrackState result;
  result.parameters=par;
  result.errors=errorProp*err*ROOT::Math::Transpose(errorProp);
  if (dump) {
    std::cout << "result.errors" << std::endl;
    dumpMatrix(result.errors);
  }
  return result;
}

TrackState propagateHelixToR_test(TrackState& inputState, int charge, float r) {

  bool dump = false;

  float xin = inputState.parameters.At(0);
  float yin = inputState.parameters.At(1);
  float zin = inputState.parameters.At(2);
  float pxin = inputState.parameters.At(3);
  float pyin = inputState.parameters.At(4);
  float pzin = inputState.parameters.At(5);

  float pt2 = pxin*pxin+pyin*pyin;
  float pt = sqrt(pt2);
  float pt3 = pt*pt2;
  //p=0.3Br => r=p/(0.3*B)
  float k=charge*100./(-0.299792458*3.8);
  float curvature = pt*k;//in cm
  if (dump) std::cout << "curvature=" << curvature << std::endl;
  float ctgTheta=pzin/pt;

  float r0in = sqrt(xin*xin+yin*yin);

  //make a copy for now...
  SVector6 par = inputState.parameters;
  SMatrix66 err = inputState.errors;

  //test//
  //try to get track circle center position
  //float xc = xin + curvature*pyin/pt;
  //float yc = yin + curvature*pxin/pt;
  //float rc = sqrt(xc*xc+yc*yc);
  //if (dump) std::cout << "rc=" << rc << " xc=" << xc << " yc=" << yc << std::endl;
  //test//
  //try to use formula in page 4 of http://www.phys.ufl.edu/~avery/fitting/fitting4.pdf
  //float c=1./(2.*curvature);
  //float D = rc-curvature;
  //float B=c*sqrt((r*r-D*D)/(1+2*c*D));//sin(0.0457164/2.);
  //float testx = xin + pxin*curvature*2.*B*sqrt(1.-B*B)/pt - pyin*curvature*2.*B*B/pt;
  //float testy = yin + pyin*curvature*2.*B*sqrt(1.-B*B)/pt + pxin*curvature*2.*B*B/pt;
  //if (dump) std::cout << "B=" << B <<  " testx=" << testx <<  " testy=" << testy << " 2*asinB=" << 2*asin(B) << std::endl;
  //test//
  //try to compute intersection between circles (approximate for small angles)
  float ceq = r0in*r0in - r*r;
  float beq = 2*k*(xin*pxin+yin*pyin);
  float aeq = k*k*(pt2+(pxin*yin-pyin*xin)/k);
  float xeq1 = (-beq + sqrt(beq*beq-4*aeq*ceq))/(2*aeq);
  float xeq2 = (-beq - sqrt(beq*beq-4*aeq*ceq))/(2*aeq);
  if (dump) std::cout << "xeq1=" << xeq1 << " xeq2=" << xeq2 << std::endl;
  //test//

  float totalAngPath=xeq1;

  float TD=totalAngPath*curvature;
  float TP=totalAngPath;
  float C=curvature;

  float cosTP = cos(TP);
  float sinTP = sin(TP);

  //fixme: these need to be derived!!!!!!!!
  float dTDdx = 0.;
  float dTDdy = 0.;
  float dTDdpx = 0.;
  float dTDdpy = 0.;
  //fixme

  par.At(0) = xin + k*(pxin*sinTP-pyin*(1-cosTP));
  par.At(1) = yin + k*(pyin*sinTP+pxin*(1-cosTP));
  par.At(2) = zin + TD*ctgTheta;
  
  par.At(3) = pxin*cosTP-pyin*sinTP;
  par.At(4) = pyin*cosTP+pxin*sinTP;
  par.At(5) = pzin;

  if (dump) std::cout << "TD=" << TD << " TP=" << TP << " arrived at r=" << sqrt(par.At(0)*par.At(0)+par.At(1)*par.At(1)) << std::endl;

  float dCdpx = k*pxin/pt;
  float dCdpy = k*pyin/pt;

  float dTPdx = dTDdx/C;
  float dTPdy = dTDdy/C;
  float dTPdpx = (dTDdpx*C - TD*dCdpx)/(C*C);
  float dTPdpy = (dTDdpy*C - TD*dCdpy)/(C*C);

  //par.At(0) = xin + k*(pxin*sinTP-pyin*(1-cosTP));
  //par.At(1) = yin + k*(pyin*sinTP+pxin*(1-cosTP));
  //par.At(2) = zin + TD*ctgTheta;

  float dxdx = 1 + k*dTPdx*(pxin*sinTP + pyin*cosTP);
  float dxdy = k*dTPdy*(pxin*sinTP + pyin*cosTP);
  float dydx = k*dTPdx*(pyin*sinTP - pxin*cosTP);
  float dydy = 1 + k*dTPdy*(pyin*sinTP - pxin*cosTP);

  float dxdpx = k*(sinTP + pxin*cosTP*dTPdpx - pyin*sinTP*dTPdpx);
  float dxdpy = k*(pxin*cosTP*dTPdpy - 1. + cosTP - pyin*sinTP*dTPdpy);
  float dydpx = k*(pyin*cosTP*dTPdpx + 1. - cosTP + pxin*sinTP*dTPdpx);
  float dydpy = k*(sinTP + pyin*cosTP*dTPdpy + pxin*sinTP*dTPdpy);

  float dzdx = dTDdx*ctgTheta;
  float dzdy = dTDdy*ctgTheta;

  float dzdpx = dTDdpx*ctgTheta - TD*pzin*pxin/pt3;
  float dzdpy = dTDdpy*ctgTheta - TD*pzin*pyin/pt3;
  float dzdpz = TD/pt;//fixme if I set this term to 0 then it works...

  //par.At(3) = pxin*cosTP-pyin*sinTP;
  //par.At(4) = pyin*cosTP+pxin*sinTP;
  //par.At(5) = pzin;

  float dpxdx = -dTPdx*(pxin*sinTP + pyin*cosTP);
  float dpxdy = -dTPdy*(pxin*sinTP + pyin*cosTP);
  float dpydx = -dTPdx*(pyin*sinTP - pxin*cosTP);
  float dpydy = -dTPdy*(pyin*sinTP - pxin*cosTP);
  
  float dpxdpx = cosTP - dTPdpx*(pxin*sinTP + pyin*cosTP);
  float dpxdpy = -sinTP - dTPdpy*(pxin*sinTP + pyin*cosTP);
  float dpydpx = +sinTP - dTPdpx*(pyin*sinTP - pxin*cosTP);
  float dpydpy = cosTP - dTPdpy*(pyin*sinTP - pxin*cosTP);

  //jacobian
  SMatrix66 errorProp = ROOT::Math::SMatrixIdentity();
  errorProp(0,0)=dxdx;
  errorProp(0,1)=dxdy;
  errorProp(0,3)=dxdpx;
  errorProp(0,4)=dxdpy;

  errorProp(1,0)=dydx;
  errorProp(1,1)=dydy;
  errorProp(1,3)=dydpx;
  errorProp(1,4)=dydpy;

  errorProp(2,0)=dzdx;
  errorProp(2,1)=dzdy;
  errorProp(2,3)=dzdpx;
  errorProp(2,4)=dzdpy;
  errorProp(2,5)=dzdpz;

  errorProp(3,0)=dpxdx;
  errorProp(3,1)=dpxdy;
  errorProp(3,3)=dpxdpx;
  errorProp(3,4)=dpxdpy;

  errorProp(4,0)=dpydx;
  errorProp(4,1)=dpydy;
  errorProp(4,3)=dpydpx;
  errorProp(4,4)=dpydpy;

  if (dump) {
    std::cout << "errorProp" << std::endl;
    dumpMatrix(errorProp);
  }

  TrackState result;
  result.parameters=par;
  result.errors=errorProp*err*ROOT::Math::Transpose(errorProp);
  if (dump) {
    std::cout << "result.errors" << std::endl;
    dumpMatrix(result.errors);
  }
  return result;

}

#endif
