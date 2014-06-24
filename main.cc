/*
  g++ -std=c++11 -O3 -Wall -Wno-unknown-pragmas -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
  icc -std=gnu++0x -O3 -openmp -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
*/

#include "fittest.h"
#include "buildtest.h"

int main(){
  bool saveTree = true;
  //runFittingTest(saveTree,5000); 
  //runFittingTestPlex(saveTree); 
  runBuildingTest(saveTree,10); 
  return 0;

}

