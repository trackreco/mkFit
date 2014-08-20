/*
  g++ -std=c++11 -O3 -Wall -Wno-unknown-pragmas -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
  icc -std=gnu++0x -O3 -openmp -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
*/

#include <sstream>

#include "Matrix.h"
#include "fittest.h"
#include "buildtest.h"
#include "Geometry.h"
#include "USolids/include/UTubs.hh"

int main(){
  Geometry* theGeom = new Geometry;

  // NB: we currently assume that each node is a layer, and that layers
  // are added starting from the center
  for (int l = 0; l < 10; l++) {
    float r = (l+1)*4;
    std::string s = "Cylinder" + std::string(1, 48+l);
    UTubs* utub = new UTubs(s, r, r+.01, 100.0, 0, TMath::TwoPi());
    theGeom->AddLayer(utub);
  }

  bool saveTree = true;
  //runFittingTest(saveTree,5000,theGeom); 
  //runFittingTestPlex(saveTree,theGeom); 
  runBuildingTest(saveTree,1,theGeom); 

  delete theGeom;
  return 0;
}
