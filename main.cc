/*
  g++ -std=c++11 -O3 -Wall -Wno-unknown-pragmas -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc ConformalUtils.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
  icc -std=gnu++0x -O3 -openmp -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
*/

#include <sstream>
#include <chrono>

#include "Matrix.h"
#include "Event.h"
#include "RootValidation.h"

#ifdef TBB
#include "tbb/task_scheduler_init.h"
#endif

//#define CYLINDER

#ifdef WITH_USOLIDS
#ifdef CYLINDER
#include "USolids/include/UTubs.hh"
#else
#include "USolids/include/UPolyhedra.hh"
#include "USolids/include/USphere.hh"
#endif

void initGeom(Geometry& geom)
{
#ifdef CYLINDER
  std::cout << "Constructing USolids Cylinder geometry" << std::endl;
#else
  std::cout << "Constructing USolids Polyhedral geometry" << std::endl;
#endif

  // NB: we currently assume that each node is a layer, and that layers
  // are added starting from the center
  float eta = 2.0; // can tune this to whatever geometry required (one can make this layer dependent as well)
  for (int l = 0; l < 10; l++) {
    float r = (l+1)*4.;
    float z = r / std::tan(2.0*std::atan(std::exp(-eta))); // calculate z extent based on eta, r
#ifdef CYLINDER
    std::string s = "Cylinder" + std::string(1, 48+l);
    UTubs* utub = new UTubs(s, r, r+.01, 100.0, 0, TMath::TwoPi());
    geom.AddLayer(utub,r,z);
#else
    float xs = 5.0; // approximate sensor size in cm
    if ( l >= 5 ) // bigger sensors in outer layers
      xs *= 2.;
    int nsectors = int(TMath::TwoPi()/(2*atan2(xs/2,r))); // keep ~constant sensors size
    std::cout << "l = " << l << ", nsectors = "<< nsectors << std::endl;
    std::string s = "PolyHedra" + std::string(1, 48+l);
    const double zPlane[] = {-z,z};
    const double rInner[] = {r,r};
    const double rOuter[] = {r+.01,r+.01};
    UPolyhedra* upolyh = new UPolyhedra(s, 0, TMath::TwoPi(), nsectors, 2, zPlane, rInner, rOuter);
    geom.AddLayer(upolyh, r, z);
#endif
  }
}
#else
void initGeom(Geometry& geom)
{
  std::cout << "Constructing SimpleGeometry Cylinder geometry" << std::endl;
  float eta = 2.0; // can tune this to whatever geometry required (one can make this layer dependent as well)

  // NB: we currently assume that each node is a layer, and that layers
  // are added starting from the center
  // NB: z is just a dummy variable, VUSolid is actually infinite in size.  *** Therefore, set it to the eta of simulation ***
  float eta = 2.0; // can tune this to whatever geometry required (one can make this layer dependent as well)
  for (int l = 0; l < 10; l++) {
    float r = (l+1)*4.;
    float z = r / std::tan(2.0*std::atan(std::exp(-eta))); // calculate z extent based on eta, r
    VUSolid* utub = new VUSolid(r, r+.01);
    geom.AddLayer(utub, r, z);
  }
}
#endif


typedef std::chrono::time_point<std::chrono::system_clock> timepoint;
typedef std::chrono::duration<double> tick;

static timepoint now()
{
  return std::chrono::system_clock::now();
}

static tick delta(timepoint& t0)
{
  timepoint t1(now());
  tick d = t1 - t0;
  t0 = t1;
  return d;
}

int main(int argc, char** argv)
{
  Geometry geom;
  initGeom(geom);
#if defined(NO_ROOT)
  Validation val;
#else
  RootValidation val("valtree.root");
#endif

  for ( unsigned int i = 0; i < geom.CountLayers(); ++i ) {
    std::cout << "Layer = " << i << ", Radius = " << geom.Radius(i) << std::endl;
  }

  unsigned int Ntracks = 500;
  unsigned int Nevents = 100;

  std::vector<tick> ticks(5);

#ifdef TBB
  auto nThread(tbb::task_scheduler_init::default_num_threads());
  if (argc > 1) {
    nThread = ::atoi(argv[1]);
  }
  std::cout << "Initializing with " << nThread << " threads." << std::endl;
  tbb::task_scheduler_init tasks(nThread);
#else
  auto nThread = 1;
#endif

  for (unsigned int evt=0; evt<Nevents; ++evt) {
    //std::cout << "EVENT #"<< evt << std::endl;
    Event ev(geom, val, nThread);

    timepoint t0(now());
    ev.Simulate(Ntracks); ticks[0] += delta(t0);
#ifdef ENDTOEND
    ev.Segment();         ticks[1] += delta(t0);
    ev.Seed();            ticks[2] += delta(t0);
    ev.Find();            ticks[3] += delta(t0);
#endif
    ev.Fit();             ticks[4] += delta(t0);
  }

  std::cout << "Ticks ";
  for (auto&& tt : ticks) {
    std::cout << tt.count() << " ";
  }
  std::cout << std::endl;

  val.saveHists();
  val.deleteHists();
  return 0;
}
