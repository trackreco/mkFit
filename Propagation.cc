#include "Propagation.h"

const double tolerance = 0.001;

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
  void propagateErrors(const HelixState& in, float totalDistance, bool dump = false);

  float x, y, z, px, py, pz;
  float k, pt, pt2, pt3, r0, curvature, ctgTheta;
  float dTDdx, dTDdy, dTDdpx, dTDdpy;
  int charge;
  TrackState& state;
};

/*  APPLE: vvsincosf(&sinAP, &cosAP, &angPath, &n); */

void HelixState::updateHelix(float distance, bool updateDeriv, bool dump)
{
  const float angPath = distance/curvature;
  if (dump) std::cout << "angPath=" << angPath << std::endl;
  const float cosAP = cos(angPath);
  const float sinAP = sin(angPath);

  //helix propagation formulas
  //http://www.phys.ufl.edu/~avery/fitting/fitting4.pdf
  SVector6& par = state.parameters;
  par.At(0) = x + k*(px*sinAP-py*(1-cosAP));
  par.At(1) = y + k*(py*sinAP+px*(1-cosAP));
  par.At(2) = z + distance*ctgTheta;
  par.At(3) = px*cosAP-py*sinAP;
  par.At(4) = py*cosAP+px*sinAP;
  par.At(5) = pz;

  if (dump) {
    std::cout << "x + " << k*(px*sinAP-py*(1-cosAP)) << std::endl;
    std::cout << "y + " << k*(py*sinAP+px*(1-cosAP)) << std::endl;
    std::cout << "z + " << distance*ctgTheta << std::endl;
    std::cout <<  "px: " << px*cosAP-py*sinAP
              << " py: " << py*cosAP+px*sinAP
              << " pz: " << pz << std::endl;
  }
  
  if (updateDeriv) {

    //update derivatives on total distance for next step, where totalDistance+=r-r0
    //now r0 depends on px and py
    const float r0inv = 1./r0;
    if (dump) std::cout << "r0=" << r0 << " r0inv=" << r0inv << " pt=" << pt << std::endl;
    //update derivative on D
    const float dAPdx = -x/(r0*curvature);
    const float dAPdy = -y/(r0*curvature);
    const float dAPdpx = -angPath*px/pt2;
    const float dAPdpy = -angPath*py/pt2;

    const float dxdx = 1 + k*dAPdx*(px*sinAP + py*cosAP);
    const float dxdy = k*dAPdy*(px*sinAP + py*cosAP);
    const float dydx = k*dAPdx*(py*sinAP - px*cosAP);
    const float dydy = 1 + k*dAPdy*(py*sinAP - px*cosAP);

    const float dxdpx = k*(sinAP + px*cosAP*dAPdpx - py*sinAP*dAPdpx);
    const float dxdpy = k*(px*cosAP*dAPdpy - 1. + cosAP - py*sinAP*dAPdpy);
    const float dydpx = k*(py*cosAP*dAPdpx + 1. - cosAP + px*sinAP*dAPdpx);
    const float dydpy = k*(sinAP + py*cosAP*dAPdpy + px*sinAP*dAPdpy);

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

void HelixState::propagateErrors(const HelixState& in, float totalDistance, bool dump)
{
  const float TD=totalDistance;
  const float TP=totalDistance/curvature;
  const float C=curvature;

  SVector6& par = state.parameters;
  if (dump) std::cout << "TD=" << TD << " TP=" << TP << " arrived at r=" << sqrt(par.At(0)*par.At(0)+par.At(1)*par.At(1)) << std::endl;

  const float dCdpx = k*in.px/pt;
  const float dCdpy = k*in.py/pt;

  const float dTPdx = dTDdx/C;
  const float dTPdy = dTDdy/C;
  const float dTPdpx = (dTDdpx*C - TD*dCdpx)/(C*C);
  const float dTPdpy = (dTDdpy*C - TD*dCdpy)/(C*C);

  const float cosTP = cos(TP);
  const float sinTP = sin(TP);

  //derive these to compute jacobian
  //x = xin + k*(pxin*sinTP-pyin*(1-cosTP));
  //y = yin + k*(pyin*sinTP+pxin*(1-cosTP));
  //z = zin + TD*ctgTheta;
  //px = pxin*cosTP-pyin*sinTP;
  //py = pyin*cosTP+pxin*sinTP;
  //pz = pzin;

  //jacobian
  SMatrix66 errorProp = ROOT::Math::SMatrixIdentity();//what is not explicitly set below is 1 (0) on (off) diagonal
  errorProp(0,0) = 1 + k*dTPdx*(in.px*sinTP + in.py*cosTP);                   //dxdx;
  errorProp(0,1) = k*dTPdy*(in.px*sinTP + in.py*cosTP);                       //dxdy;
  errorProp(0,3) = k*(sinTP + in.px*cosTP*dTPdpx - in.py*sinTP*dTPdpx);       //dxdpx;
  errorProp(0,4) = k*(in.px*cosTP*dTPdpy - 1. + cosTP - in.py*sinTP*dTPdpy);  //dxdpy;

  errorProp(1,0) = k*dTPdx*(in.py*sinTP - in.px*cosTP);                       //dydx;
  errorProp(1,1) = 1 + k*dTPdy*(in.py*sinTP - in.px*cosTP);                   //dydy;
  errorProp(1,3) = k*(in.py*cosTP*dTPdpx + 1. - cosTP + in.px*sinTP*dTPdpx);  //dydpx;
  errorProp(1,4) = k*(sinTP + in.py*cosTP*dTPdpy + in.px*sinTP*dTPdpy);       //dydpy;

  errorProp(2,0) = dTDdx*ctgTheta;                                            //dzdx;
  errorProp(2,1) = dTDdy*ctgTheta;                                            //dzdy;
  errorProp(2,3) = dTDdpx*ctgTheta - TD*in.pz*in.px/pt3;                      //dzdpx;
  errorProp(2,4) = dTDdpy*ctgTheta - TD*in.pz*in.py/pt3;                      //dzdpy;
  errorProp(2,5) = TD/pt;                                                     //dzdpz;

  errorProp(3,0) = -dTPdx*(in.px*sinTP + in.py*cosTP);                        //dpxdx;
  errorProp(3,1) = -dTPdy*(in.px*sinTP + in.py*cosTP);                        //dpxdy;
  errorProp(3,3) = cosTP - dTPdpx*(in.px*sinTP + in.py*cosTP);                //dpxdpx;
  errorProp(3,4) = -sinTP - dTPdpy*(in.px*sinTP + in.py*cosTP);               //dpxdpy;

  errorProp(4,0) = -dTPdx*(in.py*sinTP - in.px*cosTP);                        //dpydx;
  errorProp(4,1) = -dTPdy*(in.py*sinTP - in.px*cosTP);                        //dpydy;
  errorProp(4,3) = +sinTP - dTPdpx*(in.py*sinTP - in.px*cosTP);               //dpydpx;
  errorProp(4,4) = +cosTP - dTPdpy*(in.py*sinTP - in.px*cosTP);               //dpydpy;

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
TrackState propagateHelixToNextSolid(TrackState& inputState, const Geometry* theGeom) {
  bool dump = false;

  const HelixState hsin(inputState);
  TrackState result(inputState);
  HelixState hsout(result);

#ifdef CHECKSTATEVALID
  if (!hsout.state.valid) {
    return hsout.state;
  }
#endif
  if (dump) std::cout << "curvature=" << hsin.curvature << std::endl;

  float totalDistance = 0;
  auto startSolid = theGeom->InsideWhat(UVector3(hsin.x,hsin.y,hsin.z));

  // have we scattered out of the solid?
  if (hsin.r0 > 1.0 && !startSolid) {
    UVector3 here(hsin.x,hsin.y,hsin.z);
    for ( unsigned int i = 0; i < theGeom->CountLayers(); ++i ) {
      auto d = theGeom->Layer(i)->SafetyFromOutside(here, true);
      if (d < tolerance) {
        startSolid = theGeom->Layer(i);
        break;
      }
    }
    if (!startSolid) {
      std::cerr << __FILE__ << ":" << __LINE__ 
                << ": not near a solid or at origin: " << inputState.parameters
                << std::endl;
      hsout.state.valid = false;
      return hsout.state;
    }
  }

  //5 iterations is a good starting point
  const unsigned int Niter = 10;
  for (unsigned int i=0;i<Niter;++i) {
    if (dump) std::cout << "propagation iteration #" << i << std::endl;
    const float distance = std::max(theGeom->SafetyFromOutside(UVector3(hsout.x,hsout.y,hsout.z),true), tolerance);
    totalDistance += distance;

    if (dump) {
      std::cout << "r0=" << hsout.r0 << " pt=" << hsout.pt << std::endl;
      std::cout << "distance=" << distance << std::endl;
    }

    const bool updateDeriv = i+1!=Niter && hsout.r0>0.;
    hsout.updateHelix(distance, updateDeriv, dump);
    hsout.setCoords(hsout.state.parameters);

    auto currentSolid = theGeom->InsideWhat(UVector3(hsout.x,hsout.y,hsout.z));
    if (dump) {
      std::cout << "Current solid = " << currentSolid << std::endl;
    }
    if (currentSolid && currentSolid != startSolid) {
      if (dump) std::cout << "Inside next solid" << std::endl;
      break;
    }

    if ( i == (Niter-1) ) {
      std::cerr << __FILE__ << ":" << __LINE__ 
                << ": failed to converge in propagateHelixToNextSolid() after " << (i+1) << " iterations, "
                << distance 
                << ", pt = " << hsout.pt
                << std::endl;
      hsout.state.valid = false;
    }
  }

  hsout.propagateErrors(hsin, totalDistance, dump);
  return hsout.state;
}

// Propagate to the next obj
// each step travels for a path length equal to the safe step between the current position and the nearest object.
TrackState propagateHelixToLayer(TrackState& inputState, unsigned int layer, const Geometry* theGeom) {
  bool dump = false;

  const VUSolid* target = theGeom->Layer(layer);

  const HelixState hsin(inputState);
  TrackState result(inputState);
  HelixState hsout(result);

#ifdef CHECKSTATEVALID
  if (!hsout.state.valid) {
    return hsout.state;
  }
#endif

  if (theGeom->InsideWhat(UVector3(hsout.x,hsout.y,hsout.z)) == target) {
    if (dump) std::cout << "Inside target" << std::endl;
    return hsout.state;
  }
  if (dump) std::cout << "curvature=" << hsin.curvature << std::endl;

  float totalDistance = 0;

  //5 iterations is a good starting point
  const unsigned int Niter = 5;
  for (unsigned int i=0;i<Niter;++i) {
    if (dump) std::cout << "propagation iteration #" << i << std::endl;

    const float distance = std::max(target->SafetyFromOutside(UVector3(hsout.x,hsout.y,hsout.z),true), tolerance);
    totalDistance += distance;

    if (dump) {
      std::cout << "r0=" << hsout.r0 << " pt=" << hsout.pt << std::endl;
      std::cout << "distance=" << distance << std::endl;
    }

    const bool updateDeriv = i+1!=Niter && hsout.r0>0.;
    hsout.updateHelix(distance, updateDeriv, dump);
    hsout.setCoords(hsout.state.parameters);

    auto currentSolid = theGeom->InsideWhat(UVector3(hsout.x,hsout.y,hsout.z));
    if (currentSolid == target) {
      if (dump) std::cout << "Inside target" << std::endl;
      break;
    }
    if ( i == (Niter-1) ) {
      std::cerr << __FILE__ << ":" << __LINE__ 
                << ": failed to converge in propagateHelixToLayer() after " << (i+1) << " iterations, "
                << distance 
                << ", pt = " << hsout.pt
                << std::endl;
      hsout.state.valid = false;
    }
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

  const HelixState hsin(inputState);
  TrackState result(inputState);
  HelixState hsout(result);

#ifdef CHECKSTATEVALID
  if (!hsout.state.valid) {
    return hsout.state;
  }
#endif

  if (dump) {
    std::cout << "attempt propagation from r=" << hsin.r0 << " to r=" << r << std::endl
              << "x=" << hsin.x << " y=" << hsin.y << " px=" << hsin.px
              << " py=" << hsin.py << " pz=" << hsin.pz << " q=" << inputState.charge << std::endl;
  }

  if (std::abs(r-hsout.r0) < tolerance) {
    if (dump) {
      std::cout << "at target radius, returning input" << std::endl;
    }
    return hsout.state;
  }

  if (dump) std::cout << "curvature=" << hsin.curvature << std::endl;

  float totalDistance = 0;

  //5 iterations is a good starting point
  const unsigned int Niter = 10;
  for (unsigned int i=0;i<Niter;++i) {

    if (dump) std::cout << "propagation iteration #" << i << std::endl;

    const float distance = r-hsout.r0;
    totalDistance+=distance;

    if (dump) {
      std::cout << "r0=" << hsout.r0 << " pt=" << hsout.pt << std::endl
                << "distance=" << distance << std::endl;
    }
 
    bool updateDeriv = i+1!=Niter && hsout.r0>0.;
    hsout.updateHelix(distance, updateDeriv, dump);
    hsout.setCoords(hsout.state.parameters);

    if (std::abs(r-hsout.r0) < tolerance) {
      if (dump) std::cout << "distance = " << r-hsout.r0 << " at iteration=" <<  i << std::endl;
      break;
    }
    if ( i == (Niter-1) && std::abs(r-hsout.r0) > tolerance) {
      std::cerr << __FILE__ << ":" << __LINE__ 
                << ": failed to converge in propagateHelixToR() after " << (i+1) << " iterations, r = "
                << r 
                << ", hsout.r = " << hsout.r0
                << std::endl;
      hsout.state.valid = false;
    }
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


// Version with fewer temporaries and Taylor expansion of sin/cos.
// This was used to compare SMatrix / Matriplex performance.

void propagateHelixToR_fewerTemps(TrackState& inputState, float r, TrackState& result)
{
   const bool dump = false;

   float xin = inputState.parameters.At(0);
   float yin = inputState.parameters.At(1);
   float pxin = inputState.parameters.At(3);
   float pyin = inputState.parameters.At(4);
   float pzin = inputState.parameters.At(5);
   float  r0in = sqrt(xin*xin+yin*yin);
   //copy into result so that can be modified and returned
   result.parameters = inputState.parameters;
   result.errors = inputState.errors;
   result.charge = inputState.charge;
   //rename so that it is short
   SVector6& par = result.parameters;
   SMatrixSym66& err = result.errors;

#ifdef DEBUG
   if (dump) std::cout << "attempt propagation from r=" << r0in << " to r=" << r << std::endl;
   if (dump) std::cout << "x=" << xin << " y=" << yin << " px=" << pxin << " py=" << pyin << " pz=" << pzin << " q=" << inputState.charge << std::endl;
   if ((r0in-r)>=0) {
      if (dump) std::cout << "target radius same or smaller than starting point, returning input" << std::endl;
      return;
   }
#endif

   float pt2 = pxin*pxin+pyin*pyin;
   float pt = sqrt(pt2);
   float ptinv = 1./pt;
   float pt2inv = ptinv*ptinv;
   //p=0.3Br => r=p/(0.3*B)
   float k=inputState.charge*100./(-0.299792458*3.8);
   float invcurvature = 1./(pt*k);//in 1./cm
   if (dump) std::cout << "curvature=" << 1./invcurvature << std::endl;
   float ctgTheta=pzin*ptinv;
   //variables to be updated at each iterations
   //derivatives initialized to value for first iteration, i.e. distance = r-r0in
   float totalDistance = 0;
   float dTDdx = r0in>0. ? -xin/r0in : 0.;
   float dTDdy = r0in>0. ? -yin/r0in : 0.;
   float dTDdpx = 0.;
   float dTDdpy = 0.;
   //temporaries used within the loop (declare here to reduce memory operations)
   float cosAP=0.;
   float sinAP=0.;
   float dAPdx = 0.;
   float dAPdy = 0.;
   float dAPdpx = 0.;
   float dAPdpy = 0.;
   // float dxdvar = 0.;
   // float dydvar = 0.;
   //5 iterations is a good starting point
   const int Niter = 5;
   for (int i = 0; i < Niter; ++i)
   {
#ifdef DEBUG
      if (dump) std::cout << "propagation iteration #" << i << std::endl;
#endif
      float x = par.At(0);
      float y = par.At(1);
      float px = par.At(3);
      float py = par.At(4);

      float r0 = sqrt(par.At(0)*par.At(0)+par.At(1)*par.At(1));
#ifdef DEBUG
      if (dump) std::cout << "r0=" << r0 << " pt=" << pt << std::endl;
      if (dump) {
         if (r==r0) {
            std::cout << "distance = 0 at iteration=" << i << std::endl;
            break;
         }
      }
#endif
      float distance = r-r0;
      totalDistance += distance;
#ifdef DEBUG
      if (dump) std::cout << "distance=" << distance << std::endl;
#endif
      float angPath = distance*invcurvature;
#ifdef DEBUG
      if (dump) std::cout << "angPath=" << angPath << std::endl;
#endif
      // cosAP=cos(angPath);
      // sinAP=sin(angPath);
      sincos4(angPath, sinAP, cosAP);

      //helix propagation formulas
      //http://www.phys.ufl.edu/~avery/fitting/fitting4.pdf
      par.At(0) = par.At(0) + k*(px*sinAP-py*(1-cosAP));
      par.At(1) = par.At(1) + k*(py*sinAP+px*(1-cosAP));
      par.At(2) = par.At(2) + (r-r0)*ctgTheta;
      par.At(3) = px*cosAP-py*sinAP;
      par.At(4) = py*cosAP+px*sinAP;
      // par.At(5) = pz; //take this out as it is redundant

      if (i+1!=Niter && r0>0.)
      {
         //update derivatives on total distance for next step, where totalDistance+=r-r0
         //now r0 depends on px and py
         r0 = 1./r0;//WARNING, now r0 is r0inv (one less temporary)
#ifdef DEBUG
         if (dump) std::cout << "r0=" << 1./r0 << " r0inv=" << r0 << " pt=" << pt << std::endl;
#endif
         //update derivative on D
         dAPdx  = -x*r0*invcurvature;
         dAPdy  = -y*r0*invcurvature;
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
      if (dump) std::cout << par.At(0) << " " << par.At(1) << " " << par.At(2) << std::endl;
      if (dump) std::cout << par.At(3) << " " << par.At(4) << " " << par.At(5) << std::endl;
#endif
   }
   float totalAngPath=totalDistance*invcurvature;
   float& TD=totalDistance;
   float& TP=totalAngPath;
   float& iC=invcurvature;
#ifdef DEBUG
   if (dump) std::cout << "TD=" << TD << " TP=" << TP << " arrived at r=" << sqrt(par.At(0)*par.At(0)+par.At(1)*par.At(1)) << std::endl;
#endif
   float dCdpx = k*pxin*ptinv;
   float dCdpy = k*pyin*ptinv;
   float dTPdx = dTDdx*iC;
   float dTPdy = dTDdy*iC;
   float dTPdpx = (dTDdpx/iC - TD*dCdpx)*iC*iC;
   float dTPdpy = (dTDdpy/iC - TD*dCdpy)*iC*iC;

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
   SMatrix66 errorProp = ROOT::Math::SMatrixIdentity();//what is not explicitly set below is 1 (0) on (off) diagonal
   errorProp(0,0) = 1 + k*dTPdx*(pxin*sinTP + pyin*cosTP);      //dxdx;
   errorProp(0,1) = k*dTPdy*(pxin*sinTP + pyin*cosTP);  //dxdy;
   errorProp(0,3) = k*(sinTP + pxin*cosTP*dTPdpx - pyin*sinTP*dTPdpx); //dxdpx;
   errorProp(0,4) = k*(pxin*cosTP*dTPdpy - 1. + cosTP - pyin*sinTP*dTPdpy);//dxdpy;
   errorProp(1,0) = k*dTPdx*(pyin*sinTP - pxin*cosTP);  //dydx;
   errorProp(1,1) = 1 + k*dTPdy*(pyin*sinTP - pxin*cosTP);      //dydy;
   errorProp(1,3) = k*(pyin*cosTP*dTPdpx + 1. - cosTP + pxin*sinTP*dTPdpx);//dydpx;
   errorProp(1,4) = k*(sinTP + pyin*cosTP*dTPdpy + pxin*sinTP*dTPdpy); //dydpy;
   errorProp(2,0) = dTDdx*ctgTheta;     //dzdx;
   errorProp(2,1) = dTDdy*ctgTheta;     //dzdy;
   errorProp(2,3) = dTDdpx*ctgTheta - TD*pzin*pxin*pt2inv*ptinv;//dzdpx;
   errorProp(2,4) = dTDdpy*ctgTheta - TD*pzin*pyin*pt2inv*ptinv;//dzdpy;
   errorProp(2,5) = TD*ptinv; //dzdpz;
   errorProp(3,0) = -dTPdx*(pxin*sinTP + pyin*cosTP);   //dpxdx;
   errorProp(3,1) = -dTPdy*(pxin*sinTP + pyin*cosTP);   //dpxdy;
   errorProp(3,3) = cosTP - dTPdpx*(pxin*sinTP + pyin*cosTP); //dpxdpx;
   errorProp(3,4) = -sinTP - dTPdpy*(pxin*sinTP + pyin*cosTP);//dpxdpy;
   errorProp(4,0) = -dTPdx*(pyin*sinTP - pxin*cosTP); //dpydx;
   errorProp(4,1) = -dTPdy*(pyin*sinTP - pxin*cosTP);   //dpydy;
   errorProp(4,3) = +sinTP - dTPdpx*(pyin*sinTP - pxin*cosTP);//dpydpx;
   errorProp(4,4) = +cosTP - dTPdpy*(pyin*sinTP - pxin*cosTP);//dpydpy;
   result.errors=ROOT::Math::Similarity(errorProp,err);
#ifdef DEBUG
   if (dump) {
      std::cout << "errorProp" << std::endl;
      dumpMatrix(errorProp);
      std::cout << "result.errors" << std::endl;
      dumpMatrix(result.errors);
   }
#endif
   /*
     if (fabs(sqrt(par[0]*par[0]+par[1]*par[1])-r)>0.0001) {
     std::cout << "DID NOT GET TO R, dR=" << fabs(sqrt(par[0]*par[0]+par[1]*par[1])-r)
     << " r=" << r << " r0in=" << r0in << " rout=" << sqrt(par[0]*par[0]+par[1]*par[1]) << std::endl;
     std::cout << "pt=" << pt << " pz=" << inputState.parameters.At(2) << std::endl;
     }
   */
   return;
}
