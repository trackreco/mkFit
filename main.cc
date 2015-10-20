/*
  g++ -std=c++11 -O3 -Wall -Wno-unknown-pragmas -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc ConformalUtils.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
  icc -std=gnu++0x -O3 -openmp -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
*/

#include <sstream>
#include <chrono>

#include "Matrix.h"
#include "Event.h"
#include "TTreeValidation.h"

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
  for (unsigned int l = 0; l < Config::nLayers; l++) {
    float r = (l+1)*Config::fRadialSpacing;
    float z = (Config::fRadialSpacing * Config::nLayers) / std::tan(2.0*std::atan(std::exp(-Config::fEtaDet))); // calculate z extent based on eta for the last radius (used to be based on layer radius)
#ifdef CYLINDER
    std::string s = "Cylinder" + std::string(1, 48+l);
    UTubs* utub = new UTubs(s, r, r+Config::fRadialExtent, z, 0, Config::TwoPI);
    geom.AddLayer(utub,r,z);
#else
    float xs = 0.;
    if ( l < 5 ) {
      xs = Config::fInnerSensorSize;
    }
    else if ( l >= 5 ) { // bigger sensors in outer layers
      xs = Config::fOuterSensorSize;
    }
    int nsectors = int(Config::TwoPI/(2*atan2(xs/2,r))); // keep ~constant sensors size
    std::cout << "l = " << l << ", nsectors = "<< nsectors << std::endl;
    std::string s = "PolyHedra" + std::string(1, 48+l);
    const double zPlane[] = {-z,z};
    const double rInner[] = {r,r};
    const double rOuter[] = {r+Config::fRadialExtent,r+Config::fRadialExtent};
    UPolyhedra* upolyh = new UPolyhedra(s, 0, Config::TwoPI, nsectors, 2, zPlane, rInner, rOuter);
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

  for (unsigned int l = 0; l < Config::nLayers; l++) {
    float r = (l+1)*Config::fRadialSpacing;
    float z = r / std::tan(2.0*std::atan(std::exp(-Config::fEtaDet))); // calculate z extent based on eta, r
    VUSolid* utub = new VUSolid(r, r+Config::fRadialExtent);
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
  TTreeValidation val("valtree.root");
#endif

  for ( unsigned int i = 0; i < Config::nLayers; ++i ) {
    std::cout << "Layer = " << i << ", Radius = " << geom.Radius(i) << std::endl;
  }

  std::vector<tick> ticks(6);
  std::vector<unsigned int> tracks(4);
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

  for (unsigned int evt=0; evt<Config::nEvents; ++evt) {
    Event ev(geom, val, evt, nThread);
    std::cout << "EVENT #"<< ev.evtID() << std::endl;

    timepoint t0(now());
#ifdef ENDTOEND
    ev.Simulate();           ticks[0] += delta(t0);
    ev.Segment();            ticks[1] += delta(t0);
    ev.Seed();               ticks[2] += delta(t0);
    ev.Find();               ticks[3] += delta(t0);
    ev.Fit();                ticks[4] += delta(t0);
    ev.Validate(ev.evtID()); ticks[5] += delta(t0);
#endif
    std::cout << "sim: " << ev.simTracks_.size() << " seed: " << ev.seedTracks_.size() << " found: " << ev.candidateTracks_.size() << " fit: " << ev.fitTracks_.size() << std::endl;
    tracks[0] += ev.simTracks_.size();
    tracks[1] += ev.seedTracks_.size();
    tracks[2] += ev.candidateTracks_.size();
    tracks[3] += ev.fitTracks_.size();
  }

  std::vector<double> time(ticks.size());
  for (unsigned int i = 0; i < ticks.size(); ++i){
    time[i]=ticks[i].count();
  }

  val.fillConfigTree(time);
  val.saveTTrees(); 

  std::cout << "Ticks ";
  for (auto&& tt : time) {
    std::cout << tt << " ";
  }
  std::cout << std::endl;
  std::cout << "Tracks ";
  for (auto t : tracks) {
    std::cout << t << " ";
  }
  std::cout << std::endl;
  return 0;
}
