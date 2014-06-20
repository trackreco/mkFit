#include "ConformalUtils.h"

//M. Hansroul, H. Jeremie and D. Savard, NIM A 270 (1988) 498
//http://www.sciencedirect.com/science/article/pii/016890028890722X
void conformalFit(Hit hit0,Hit hit1,Hit hit2,int charge) {

  float bsx = 0.;
  float bsy = 0.;

  float x[3],y[3];
  x[0]=hit0.position()[0]-bsx;
  x[1]=hit1.position()[0]-bsx;
  x[2]=hit2.position()[0]-bsx;
  y[0]=hit0.position()[1]-bsy;
  y[1]=hit1.position()[1]-bsy;
  y[2]=hit2.position()[1]-bsy;

  float u[3],v[3];
  for (unsigned int h=0;h<3;++h) {
    u[h]=x[h]/(x[h]*x[h]+y[h]*y[h]);
    v[h]=y[h]/(x[h]*x[h]+y[h]*y[h]);
  }

  // R^2=a^2+b^2
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
  float R=sqrt(a*a+b*b);
  float e=b*b*b*C[2]/(R*R*R);

  float k=charge*100./(-0.299792458*3.8);
  float pt = R/k;
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

}
