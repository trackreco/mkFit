#include "ConformalUtils.h"

//M. Hansroul, H. Jeremie and D. Savard, NIM A 270 (1988) 498
//http://www.sciencedirect.com/science/article/pii/016890028890722X

void conformalFit(const Hit& hit0, const Hit& hit1, const Hit& hit2, int charge, TrackState& fitStateHit0, bool fiterrs) {

  // store hit info
  float x[3],y[3],z[3];
  x[0]=hit0.position()[0];
  x[1]=hit1.position()[0];
  x[2]=hit2.position()[0];
  y[0]=hit0.position()[1];
  y[1]=hit1.position()[1];
  y[2]=hit2.position()[1];
  z[0]=hit0.position()[2];
  z[1]=hit1.position()[2];
  z[2]=hit2.position()[2];

  // save x,y,z from first hit
  fitStateHit0.parameters[0] = x[0];
  fitStateHit0.parameters[1] = y[0];
  fitStateHit0.parameters[2] = z[0];

  //r-phi smearing
  float hitRad2   = x[0]*x[0]+y[0]*y[0];
  float varPosPhi = Config::varXY/hitRad2;

  //this gives nice fitting results when scaling the errors by 10
  fitStateHit0.errors=ROOT::Math::SMatrixIdentity();

  //r-phi smearing
  fitStateHit0.errors[0][0] = x[0]*x[0]*Config::varR/hitRad2 + y[0]*y[0]*varPosPhi;
  fitStateHit0.errors[1][1] = x[0]*x[0]*varPosPhi + y[0]*y[0]*Config::varR/hitRad2;
  fitStateHit0.errors[2][2] = Config::varZ;
  fitStateHit0.errors[1][0] = x[0]*y[0]*(Config::varR/hitRad2 - varPosPhi);
  fitStateHit0.errors[0][1] = fitStateHit0.errors[1][0];

  // now do calculation of momenta
  float u[3],v[3]; // conformal points
  float initphi = fabs(getPhi(x[1],y[1])); // use this to decide when to use x -> u or x -> v
  bool xtou = (initphi<Config::PIOver4 || initphi>Config::PI3Over4);

  if (xtou){ // x -> u 
    for (unsigned int h=0;h<3;++h) {
      u[h]=x[h]/getRad2(x[h],y[h]);
      v[h]=y[h]/getRad2(x[h],y[h]);
    }
  }
  else { // x -> v
    for (unsigned int h=0;h<3;++h) {
      v[h]=x[h]/getRad2(x[h],y[h]);
      u[h]=y[h]/getRad2(x[h],y[h]);
    }
  }

  //R^2=a^2+b^2
  //v = 1/2b - ua/b - u^2*e*R^3/b^3
  SVector3 B(v[0],v[1],v[2]);
  SMatrix33 A;
  for (unsigned int h=0;h<3;++h) {
    A(h,0) = 1.;A(h,1) = -u[h]; A(h,2) = -u[h]*u[h];
  }
  A.Invert();
  SVector3 C=A*B;

  const float b=1./(2.*C[0]);
  const float a=b*C[1];

  // Special note, "vr" is short for vector, not vertex!           
  // evaluate momentum, phi, theta at layer one
  // taking vector from center of circle to layer one position
  // therefore phi is the perpendicular to vector just described        

  const float vrx = (xtou ? x[0]-a : x[0]-b);
  const float vry = (xtou ? y[0]-b : y[0]-a);

  const float R   = sqrtf(getRad2(vrx,vry));
  //float e=b*b*b*C[2]/(R*R*R);
  const float k   = charge*100./(-Config::sol*Config::Bfield);
  // compute pt
  const float pt  = R/k;
  // compute phi
  const float phi = getPhi(vry,vrx);
  //compute theta
  const float tantheta = sqrtf(getRad2((x[0]-x[2]),(y[0]-y[2])))/(z[2]-z[0]);

#ifdef CCSCOORD
  fitStateHit0.parameters[3] = 1.0f/pt;
  fitStateHit0.parameters[4] = phi;
  fitStateHit0.parameters[5] = std::tan(tantheta);
#ifdef INWARDFIT
  if (fiterrs) fitStateHit0.parameters[5] *= -1.0f; // tangent is an odd function, so tan(pt/pz) when inward is tan(pt/-pz)= -tan(pt/pz)
#endif
  fitStateHit0.errors[3][3] = (fiterrs ? Config::ptinverr049 * Config::ptinverr049 : Config::ptinverr012 * Config::ptinverr012);
  fitStateHit0.errors[4][4] = (fiterrs ? Config::phierr049   * Config::phierr049   : Config::phierr012   * Config::phierr012);
  fitStateHit0.errors[5][5] = (fiterrs ? Config::thetaerr049 * Config::thetaerr049 : Config::thetaerr012 * Config::thetaerr012);
#else
  float px = fabs(pt*cos(phi))*((x[1]-x[0])>0. ? 1. : -1.);
  float py = fabs(pt*sin(phi))*((y[1]-y[0])>0. ? 1. : -1.);
  float pz = fabs(pt/tantheta)*((z[1]-z[0])>0. ? 1. : -1.);
#ifdef INWARDFIT
  if (fiterrs) { // need conformal fit on seeds to be forward!
    px*=-1.;
    py*=-1.;
    pz*=-1.;
  }
#endif
  //return px,py,pz
  fitStateHit0.parameters[3] = px;
  fitStateHit0.parameters[4] = py;
  fitStateHit0.parameters[5] = pz;

  //get them a posteriori from width of residue plots (i.e. unitary pulls) + global maxima scan of nHits / track in super debug mode
  const float pt2 = pt*pt;
  const float pz2 = pz*pz;

  const float varPt    = (fiterrs ? Config::ptinverr049 * Config::ptinverr049 : Config::ptinverr012 * Config::ptinverr012) * pt2;
  const float varPhi   = (fiterrs ? Config::phierr049   * Config::phierr049   : Config::phierr012   * Config::phierr012);
  const float varTheta = (fiterrs ? Config::thetaerr049 * Config::thetaerr049 : Config::thetaerr012 * Config::thetaerr012);

  fitStateHit0.errors[3][3] = px*px*varPt + py*py*varPhi;
  fitStateHit0.errors[4][4] = py*py*varPt + px*px*varPhi;
  fitStateHit0.errors[5][5] = pz2  *varPt + ((pz2+pt2)*(pz2+pz2)/pt2)*varTheta;
#endif // cartesian coords

  fitStateHit0.charge = charge; //taken from slopes!
  //dumpMatrix(fitStateHit0.errors);
}
