#include "ConformalUtils.h"

//M. Hansroul, H. Jeremie and D. Savard, NIM A 270 (1988) 498
//http://www.sciencedirect.com/science/article/pii/016890028890722X

void conformalFit(const Hit& hit0, const Hit& hit1, const Hit& hit2, int charge, TrackState& fitStateHit0, bool fiterrs) {

  //fixme: does this work in case bs not in (0,0)? I think so, but need to check

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

  float u[3],v[3];
  for (unsigned int h=0;h<3;++h) {
    u[h]=x[h]/(x[h]*x[h]+y[h]*y[h]);
    v[h]=y[h]/(x[h]*x[h]+y[h]*y[h]);
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

  float b=1./(2.*C[0]);
  float a=b*C[1];
  float R=sqrt((x[0]-a)*(x[0]-a)+(y[0]-b)*(y[0]-b));
  //float e=b*b*b*C[2]/(R*R*R);

  float k=charge*100./(-Config::sol*Config::Bfield);
  float pt = R/k;
  /*
  std::cout << "hit0=" << x[0] << "," << y[0] << std::endl;
  std::cout << "hit1=" << x[1] << "," << y[1] << std::endl;
  std::cout << "hit2=" << x[2] << "," << y[2] << std::endl;
  std::cout << "hit0t=" << u[0] << "," << v[0] << std::endl;
  std::cout << "hit1t=" << u[1] << "," << v[1] << std::endl;
  std::cout << "hit2t=" << u[2] << "," << v[2] << std::endl;
  std::cout << "vfit0=" << C[0]-u[0]*C[1]-u[0]*u[0]*C[2] << std::endl;
  std::cout << "vfit1=" << C[0]-u[1]*C[1]-u[1]*u[1]*C[2] << std::endl;
  std::cout << "vfit2=" << C[0]-u[2]*C[1]-u[2]*u[2]*C[2] << std::endl;
  std::cout << "c0=" << C[0] << " c1=" << C[1] << " c2=" << C[2] << std::endl;
  std::cout << "a=" << a << " b=" << b << " e=" << e << std::endl;
  std::cout << "R=" << R << " pt=" << pt << std::endl;
  */

  float vrx = x[0]-a;
  float vry = y[0]-b;
  float phi = atan2(vrx,vry);
  float px = fabs(pt*cos(phi))*((x[1]-x[0])>0. ? 1. : -1.);
  float py = fabs(pt*sin(phi))*((y[1]-y[0])>0. ? 1. : -1.);

  //compute theta
  float tantheta = sqrt((x[0]-x[2])*(x[0]-x[2])+(y[0]-y[2])*(y[0]-y[2]))/(z[2]-z[0]);
  float pz = fabs(pt/tantheta)*((z[1]-z[0])>0. ? 1. : -1.);
#ifdef INWARDFIT
  if (fiterrs) { // need conformal fit on seeds to be forward!
    px*=-1.;
    py*=-1.;
    pz*=-1.;
  }
#endif
  //return px,py,pz
  //std::cout << "fit px=" << px << " py=" << py << " pz=" << pz << std::endl; 
  fitStateHit0.parameters[0] = x[0];
  fitStateHit0.parameters[1] = y[0];
  fitStateHit0.parameters[2] = z[0];
  fitStateHit0.parameters[3] = px;
  fitStateHit0.parameters[4] = py;
  fitStateHit0.parameters[5] = pz;
  //get them a posteriori from width of residue plots (i.e. unitary pulls)

  //r-phi smearing
  float hitRad2 = x[0]*x[0]+y[0]*y[0];
  float varPhi = Config::varXY/hitRad2;

  float ptinverr = 0.;
  float pterr    = 0.;
  float phierr   = 0.;
  float thetaerr = 0.;

  if (fiterrs) { // use fit errors, ie. hits on layers 0,5,9
    ptinverr = Config::ptinverr049;
    phierr   = Config::phierr049;
    thetaerr = Config::thetaerr049;
  }
  else{ //use seed errors, ie. hits on layers 0,1,2
    ptinverr = Config::ptinverr012;
    phierr   = Config::phierr012;
    thetaerr = Config::thetaerr012;
  }
  pterr = (pt*pt)*ptinverr;
   
  //this gives nice fitting results when scaling the errors by 10
  fitStateHit0.errors=ROOT::Math::SMatrixIdentity();
  //xy smearing
  //fitStateHit0.errors[0][0] = pow(xerr,2);
  //fitStateHit0.errors[1][1] = pow(yerr,2);
  //fitStateHit0.errors[2][2] = pow(zerr,2);
  //r-phi smearing
  fitStateHit0.errors[0][0] = x[0]*x[0]*Config::varR/hitRad2 + y[0]*y[0]*varPhi;
  fitStateHit0.errors[1][1] = x[0]*x[0]*varPhi + y[0]*y[0]*Config::varR/hitRad2;
  fitStateHit0.errors[2][2] = Config::varZ;
  fitStateHit0.errors[1][0] = x[0]*y[0]*(Config::varR/hitRad2 - varPhi);
  fitStateHit0.errors[0][1] = fitStateHit0.errors[1][0];

  fitStateHit0.errors[3][3] = pow(cos(phi),2)*pow(pterr,2)+pow(pt*sin(phi),2)*pow(phierr,2);
  fitStateHit0.errors[4][4] = pow(sin(phi),2)*pow(pterr,2)+pow(pt*cos(phi),2)*pow(phierr,2);
  fitStateHit0.errors[5][5] = pow(1./tantheta,2)*pow(pterr,2)+pow(pt/pow(sin(atan(tantheta)),2),2)*pow(thetaerr,2);

  /*
  //fixme: if done with correlations pt pull gets larger, do I have a bug?  (actually scaling by 10k it looks nice as well)
  SMatrixSym66 fiterrors = ROOT::Math::SMatrixIdentity();//x,y,z,pt,phi,theta
  fiterrors[0][0] = xerr*xerr;
  fiterrors[1][1] = yerr*yerr;
  fiterrors[2][2] = zerr*zerr;
  fiterrors[3][3] = pterr*pterr;
  fiterrors[4][4] = phierr*phierr;
  fiterrors[5][5] = thetaerr*thetaerr;
  float dpxdpt = cos(phi);
  float dpxdphi = -py;
  float dpydpt = sin(phi);
  float dpydphi = px;
  float dpzdpt = 1/tantheta;
  float dpzdtheta = -pt/pow(sin(atan2(pt,pz)),2);
  SMatrix66 jacobian = ROOT::Math::SMatrixIdentity();//from x,y,z,pt,phi,theta to x,y,z,px,py,pz
  jacobian[3][3] = dpxdpt;
  jacobian[3][4] = dpxdphi;
  jacobian[4][3] = dpydpt;
  jacobian[4][4] = dpydphi;
  jacobian[5][3] = dpzdpt;
  jacobian[5][5] = dpzdtheta;
  fitStateHit0.errors = ROOT::Math::Similarity(jacobian,fiterrors);
  */

  fitStateHit0.charge = charge;//fixme, estimate from fit
  //dumpMatrix(fitStateHit0.errors);
}
