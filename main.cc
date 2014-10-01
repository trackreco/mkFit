/*
  g++ -std=c++11 -O3 -Wall -Wno-unknown-pragmas -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc ConformalUtils.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
  icc -std=gnu++0x -O3 -openmp -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
*/

#include <sstream>

#include "Matrix.h"
#include "fittest.h"
#include "buildtest.h"
#include "Geometry.h"
#include "USolids/include/UTubs.hh"
#include "USolids/include/UPolyhedra.hh"
#include "USolids/include/USphere.hh"

int main(){
  Geometry* theGeom = new Geometry;

  // NB: we currently assume that each node is a layer, and that layers
  // are added starting from the center
  for (int l = 0; l < 10; l++) {
    float r = (l+1)*4.;
    //#define CYLINDER
#ifdef CYLINDER
    std::string s = "Cylinder" + std::string(1, 48+l);
    UTubs* utub = new UTubs(s, r, r+.01, 100.0, 0, TMath::TwoPi());
    theGeom->AddLayer(utub,r);
#else
    float xs = 5.0; // approximate sensor size in cm
    if ( l >= 5 ) // bigger sensors in outer layers
      xs *= 2.;
    int nsectors = int(TMath::TwoPi()/(2*atan2(xs/2,r))); // keep ~constant sensors size
    std::cout << "l = " << l << ", nsectors = "<< nsectors << std::endl;
    std::string s = "PolyHedra" + std::string(1, 48+l);
    const double zPlane[] = {-100.,100.};
    const double rInner[] = {r,r};
    const double rOuter[] = {r+.01,r+.01};
    UPolyhedra* upolyh = new UPolyhedra(s, 0, TMath::TwoPi(), nsectors, 2, zPlane, rInner, rOuter);
    theGeom->AddLayer(upolyh, r);
#endif
  }

  for ( int i = 0; i < theGeom->CountLayers(); ++i ) {
    std::cout << "Layer = " << i << ", Radius = " << theGeom->Radius(i) << std::endl;
  }

  bool saveTree = true;
  runFittingTest(saveTree,50000,theGeom); 
  //runFittingTestPlex(saveTree,theGeom); 
  runBuildingTest(saveTree,100,theGeom); 

  delete theGeom;
  return 0;
}
