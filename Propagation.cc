#include "Propagation.h"
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#endif

// line propagation from state radius to hit radius
// assuming radial direction (i.e. origin at (0,0))
TrackState propagateLineToR(TrackState& inputState, float r) {
  bool dump = false;

  SVector6& par = inputState.parameters;
  SMatrixSym66& err = inputState.errors;

  //straight line for now
  float r0 = sqrt(par.At(0)*par.At(0)+par.At(1)*par.At(1));
  float dr = r-r0;
  float pt = sqrt(par.At(3)*par.At(3)+par.At(4)*par.At(4));
  float path = dr/pt;//this works only if direction is along radius, i.e. origin is at (0,0)

  TrackState result;
  result.charge = inputState.charge;

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

  result.errors=ROOT::Math::Similarity(propMatrix,err);
  return result;
}

struct HelixState {
  HelixState(TrackState& s) : state(s) {
    setCoords(s.parameters);
    setHelixPar(s);
  }

  void setCoords(SVector6& par) {
    x = par.At(0);
    y = par.At(1);
    z = par.At(2);
    px = par.At(3);
    py = par.At(4);
    pz = par.At(5);
    r0 = sqrt(x*x+y*y);
  }

  void setHelixPar(TrackState& s) {
    charge = s.charge;

    pt2 = px*px+py*py;
    pt = sqrt(pt2);
    pt3 = pt*pt2;

    //p=0.3Br => r=p/(0.3*B)
    k = charge*100./(-0.299792458*3.8);
    curvature = pt*k; //in cm
    ctgTheta=pz/pt;

    //variables to be updated at each iterations
    //derivatives initialized to value for first iteration, i.e. distance = r-r0in
    dTDdx = r0>0. ? -x/r0 : 0.;
    dTDdy = r0>0. ? -y/r0 : 0.;
    dTDdpx = 0.;
    dTDdpy = 0.;
  }

  void updateHelix(float distance, bool updateDeriv, bool dump = false);
  void propagateErrors(HelixState& in, float totalDistance, bool dump = false);

  float x, y, z, px, py, pz;
  float k, pt, pt2, pt3, r0, curvature, ctgTheta;
  float dTDdx, dTDdy, dTDdpx, dTDdpy;
  int charge;
  TrackState& state;
};

void HelixState::updateHelix(float distance, bool updateDeriv, bool dump)
{
  float angPath = distance/curvature;
  if (dump) std::cout << "angPath=" << angPath << std::endl;
  float cosAP = cos(angPath);
  float sinAP = sin(angPath);
#ifdef __APPLE__
  int n = 1;
  vvsincosf(&sinAP, &cosAP, &angPath, &n);
#else
  cosAP = cos(angPath);
  sinAP = sin(angPath);
#endif

  //helix propagation formulas
  //http://www.phys.ufl.edu/~avery/fitting/fitting4.pdf
  SVector6& par = state.parameters;
  par.At(0) = x + k*(px*sinAP-py*(1-cosAP));
  par.At(1) = y + k*(py*sinAP+px*(1-cosAP));
  par.At(2) = z + distance*ctgTheta;
  par.At(3) = px*cosAP-py*sinAP;
  par.At(4) = py*cosAP+px*sinAP;
  par.At(5) = pz;
  
  if (updateDeriv) {

    //update derivatives on total distance for next step, where totalDistance+=r-r0
    //now r0 depends on px and py
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

  if (dump) {
    std::cout << par.At(0) << " " << par.At(1) << " " << par.At(2) << std::endl
              << par.At(3) << " " << par.At(4) << " " << par.At(5) << std::endl;
  }
}

void HelixState::propagateErrors(HelixState& in, float totalDistance, bool dump)
{
  float totalAngPath=totalDistance/curvature;

  float& TD=totalDistance;
  float& TP=totalAngPath;
  float& C=curvature;

  SVector6& par = state.parameters;
  if (dump) std::cout << "TD=" << TD << " TP=" << TP << " arrived at r=" << sqrt(par.At(0)*par.At(0)+par.At(1)*par.At(1)) << std::endl;

  float dCdpx = k*in.px/pt;
  float dCdpy = k*in.py/pt;

  float dTPdx = dTDdx/C;
  float dTPdy = dTDdy/C;
  float dTPdpx = (dTDdpx*C - TD*dCdpx)/(C*C);
  float dTPdpy = (dTDdpy*C - TD*dCdpy)/(C*C);

  float cosTP = cos(TP);
  float sinTP = sin(TP);

  //derive these to compute jacobian
  //x = xin + k*(pxin*sinTP-pyin*(1-cosTP));
  //y = yin + k*(pyin*sinTP+pxin*(1-cosTP));
  //z = zin + TD*ctgTheta;
  //px = pxin*cosTP-pyin*sinTP;
  //py = pyin*cosTP+pxin*sinTP;
  //pz = pzin;

  //jacobian
  SMatrix66 errorProp = ROOT::Math::SMatrixIdentity();//what is not explicitly set below is 1 (0) on (off) diagonal
  errorProp(0,0) = 1 + k*dTPdx*(in.px*sinTP + in.py*cosTP);	          //dxdx;
  errorProp(0,1) = k*dTPdy*(in.px*sinTP + in.py*cosTP);		          //dxdy;
  errorProp(0,3) = k*(sinTP + in.px*cosTP*dTPdpx - in.py*sinTP*dTPdpx);     //dxdpx;
  errorProp(0,4) = k*(in.px*cosTP*dTPdpy - 1. + cosTP - in.py*sinTP*dTPdpy);//dxdpy;

  errorProp(1,0) = k*dTPdx*(in.py*sinTP - in.px*cosTP);		          //dydx;
  errorProp(1,1) = 1 + k*dTPdy*(in.py*sinTP - in.px*cosTP);	          //dydy;
  errorProp(1,3) = k*(in.py*cosTP*dTPdpx + 1. - cosTP + in.px*sinTP*dTPdpx);//dydpx;
  errorProp(1,4) = k*(sinTP + in.py*cosTP*dTPdpy + in.px*sinTP*dTPdpy);     //dydpy;

  errorProp(2,0) = dTDdx*ctgTheta;		      //dzdx;
  errorProp(2,1) = dTDdy*ctgTheta;		      //dzdy;
  errorProp(2,3) = dTDdpx*ctgTheta - TD*in.pz*in.px/pt3;//dzdpx;
  errorProp(2,4) = dTDdpy*ctgTheta - TD*in.pz*in.py/pt3;//dzdpy;
  errorProp(2,5) = TD/pt;                             //dzdpz;

  errorProp(3,0) = -dTPdx*(in.px*sinTP + in.py*cosTP);	     //dpxdx;
  errorProp(3,1) = -dTPdy*(in.px*sinTP + in.py*cosTP);	     //dpxdy;
  errorProp(3,3) = cosTP - dTPdpx*(in.px*sinTP + in.py*cosTP); //dpxdpx;
  errorProp(3,4) = -sinTP - dTPdpy*(in.px*sinTP + in.py*cosTP);//dpxdpy;

  errorProp(4,0) = -dTPdx*(in.py*sinTP - in.px*cosTP);         //dpydx;
  errorProp(4,1) = -dTPdy*(in.py*sinTP - in.px*cosTP);	     //dpydy;
  errorProp(4,3) = +sinTP - dTPdpx*(in.py*sinTP - in.px*cosTP);//dpydpx;
  errorProp(4,4) = +cosTP - dTPdpy*(in.py*sinTP - in.px*cosTP);//dpydpy;

  state.errors=ROOT::Math::Similarity(errorProp,state.errors);

  if (dump) {
    std::cout << "errorProp" << std::endl;
    dumpMatrix(errorProp);
    std::cout << "result.errors" << std::endl;
    dumpMatrix(state.errors);
  }
}

// helix propagation in steps along helix trajectory, several versions
// for track with pT>=1 GeV this converges to the correct path lenght in <5 iterations
// derivatives need to be updated at each iteration

// Propagate to the next obj
// each step travels for a path length equal to the safe step between the current position and the nearest object.
TrackState propagateHelixToNextSolid(TrackState& inputState, Geometry* theGeom) {
  bool dump = false;

  HelixState hsin(inputState);
  HelixState hsout(hsin);

  if (dump) std::cout << "curvature=" << hsin.curvature << std::endl;

  float totalDistance = 0;
  auto startSolid = theGeom->InsideWhat(UVector3(hsin.x,hsin.y,hsin.z));

  if (dump) {
    std::cout << "startSolid = " << startSolid << std::endl;
    // distance to all layers
    UVector3 origin(0,0,0);
    //UVector3 direction(par(4),par(5),par(6));
    UVector3 direction(hsout.px,hsout.py,hsout.pz);
    direction.Normalize();
    std::cout << "direction is " << direction << std::endl;
    for ( unsigned int i = 0; i < theGeom->CountLayers(); ++i ) {
      auto distance1 = theGeom->Layer(i)->DistanceToIn(origin, direction);
      auto distance2 = theGeom->Layer(i)->SafetyFromOutside(origin, false);
      std::cout << "disance to layer " << i << " = " << distance1 
		<< " " << distance2 << ", inside=" << theGeom->Layer(i)->Inside(origin)
		<< std::endl;
    }
  }
    

  //5 iterations is a good starting point
  unsigned int Niter = 10;
  for (unsigned int i=0;i<Niter;++i) {
    if (dump) std::cout << "propagation iteration #" << i << std::endl;
    hsout.setCoords(hsout.state.parameters);

    auto currentSolid = theGeom->InsideWhat(UVector3(hsout.x,hsout.y,hsout.z));
    if ( dump ) 
      std::cout << "Current solid = " << currentSolid << std::endl;
    if (currentSolid && currentSolid != startSolid) {
      if (dump) std::cout << "Inside next solid" << std::endl;
      break;
    }

    float distance = std::max(theGeom->SafetyFromOutside(UVector3(hsout.x,hsout.y,hsout.z),true), .0001);
    totalDistance += distance;

    if (dump) {
      std::cout << "r0=" << hsout.r0 << " pt=" << hsout.pt << std::endl;
      std::cout << "distance=" << distance << std::endl;
    }

    bool updateDeriv = i+1!=Niter && hsout.r0>0.;
    hsout.updateHelix(distance, updateDeriv, dump);
  }

  hsout.propagateErrors(hsin, totalDistance, dump);
  return hsout.state;
}

// Propagate to the next obj
// each step travels for a path length equal to the safe step between the current position and the nearest object.
TrackState propagateHelixToLayer(TrackState& inputState, unsigned int layer, Geometry* theGeom) {
  bool dump = false;

  const VUSolid* target = theGeom->Layer(layer);

  HelixState hsin(inputState);
  HelixState hsout(hsin);

  if (dump) std::cout << "curvature=" << hsin.curvature << std::endl;

  float totalDistance = 0;

  //5 iterations is a good starting point
  unsigned int Niter = 5;
  for (unsigned int i=0;i<Niter;++i) {
    if (dump) std::cout << "propagation iteration #" << i << std::endl;
    hsout.setCoords(hsout.state.parameters);

    auto currentSolid = theGeom->InsideWhat(UVector3(hsout.x,hsout.y,hsout.z));
    if (currentSolid == target) {
      if (dump) std::cout << "Inside target" << std::endl;
      break;
    }

    float distance = std::max(target->SafetyFromOutside(UVector3(hsout.x,hsout.y,hsout.z),true), .0001);
    totalDistance += distance;

    if (dump) {
      std::cout << "r0=" << hsout.r0 << " pt=" << hsout.pt << std::endl;
      std::cout << "distance=" << distance << std::endl;
    }

    bool updateDeriv = i+1!=Niter && hsout.r0>0.;
    hsout.updateHelix(distance, updateDeriv, dump);
  }

  hsout.propagateErrors(hsin, totalDistance, dump);
  return hsout.state;
}


// helix propagation in steps along helix trajectory. 
// each step travels for a path lenght equal to delta r between the current position and the target radius. 
// for track with pT>=1 GeV this converges to the correct path lenght in <5 iterations
// derivatives need to be updated at each iteration
TrackState propagateHelixToR(TrackState& inputState, float r) {
  bool dump = false;

  HelixState hsin(inputState);
  HelixState hsout(hsin);

  if (dump) {
    std::cout << "attempt propagation from r=" << hsin.r0 << " to r=" << r << std::endl
              << "x=" << hsin.x << " y=" << hsin.y << " px=" << hsin.px
              << " py=" << hsin.py << " pz=" << hsin.pz << " q=" << inputState.charge << std::endl;
  }

  if ((hsin.r0-r)>=0) {
    if (dump) std::cout << "target radius same or smaller than starting point, returning input" << std::endl;
    return hsin.state;
  }

  if (dump) std::cout << "curvature=" << hsin.curvature << std::endl;

  float totalDistance = 0;

  //5 iterations is a good starting point
  unsigned int Niter = 5;
  for (unsigned int i=0;i<Niter;++i) {

    if (dump) std::cout << "propagation iteration #" << i << std::endl;
    hsout.setCoords(hsout.state.parameters);

    if (r==hsout.r0) {
      if (dump) std::cout << "distance = 0 at iteration=" <<  i << std::endl;
      break;
    }

    float distance = r-hsout.r0;
    totalDistance+=distance;

    if (dump) {
      std::cout << "r0=" << hsout.r0 << " pt=" << hsout.pt << std::endl
                << "distance=" << distance << std::endl;
    }
 
    bool updateDeriv = i+1!=Niter && hsout.r0>0.;
    hsout.updateHelix(distance, updateDeriv, dump);
  }

  hsout.propagateErrors(hsin, totalDistance, dump);
  return hsout.state;
}

//test towards a helix propagation without iterative approach
//version below solves the equation for the angular path at which x^2+y^2=r^2
//problems: 1. need first order approximation of sin and cos, 
//2. there are 2 numerical solutions, 3. need to propagate uncertainties throgh the 2nd order equation
TrackState propagateHelixToR_test(TrackState& inputState, float r) {

  bool dump = false;

  int charge = inputState.charge;

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

  //make a copy so that can be modified and returned
  SVector6 par = inputState.parameters;
  SMatrixSym66 err = inputState.errors;

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
  //solve 2nd order equation, obtained setting x^2+y^2=r^2 and solve for the angular path
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
  result.errors=ROOT::Math::Similarity(errorProp,err);
  result.charge = charge;
  if (dump) {
    std::cout << "result.errors" << std::endl;
    dumpMatrix(result.errors);
  }
  return result;

}
