/*
  g++ -std=c++11 -O3 -Wall -Wno-unknown-pragmas -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc ConformalUtils.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
  icc -std=gnu++0x -O3 -openmp -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
*/

#include <sstream>

#include "Matrix.h"
#include "Event.h"
#include "RootValidation.h"

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

  // NB: we currently assume that each node is a layer, and that layers
  // are added starting from the center
  // NB: z is just a dummy variable, VUSolid is actually infinite in size.  *** Therefore, set it to the eta of simulation ***
  for (int l = 0; l < 10; l++) {
    float r = (l+1)*4.;
    VUSolid* utub = new VUSolid(r, r+.01);
    geom.AddLayer(utub, r, z);
  }
}
#endif


typedef long long tick_t;

static tick_t now()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return static_cast<long long>(1000000)*static_cast<long long>(tv.tv_sec) + static_cast<long long>(tv.tv_usec);
}

static tick_t delta(tick_t& t0)
{
  tick_t t1 = now();
  tick_t d = t1 - t0;
  t0 = t1;
  return d;
}

int main(){
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

  tick_t tick[5] = {0,0,0,0,0};

  for (unsigned int evt=0; evt<Nevents; ++evt) {
    std::cout << "EVENT #"<< evt << std::endl;
    Event ev(geom, val);

    tick_t t0(now());
    ev.Simulate(Ntracks); tick[0] += delta(t0);
    ev.Segment();         tick[1] += delta(t0);
    ev.Seed();            tick[2] += delta(t0);
    ev.Find();            tick[3] += delta(t0);
    ev.Fit();             tick[4] += delta(t0);
  }

  std::cout << "Ticks ";
  for (auto&& tt : tick) {
    std::cout << tt/1000000.0 << " ";
  }
  std::cout << std::endl;

  val.saveHists();
  val.deleteHists();
  return 0;
}
