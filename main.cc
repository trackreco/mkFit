/*
  g++ -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
  icc -std=gnu++0x -O3 -openmp -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
*/

#include <iostream>

#include "fittest.h"
#include "buildtest.h"

int main() {

  bool saveTree = true;

  //runFittingTest(saveTree);
  runFittingTestPlex(saveTree);
  //runBuildingTest(saveTree,10);
  return 0;

}

