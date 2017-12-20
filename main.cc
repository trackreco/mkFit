/*
  g++ -std=c++11 -O3 -Wall -Wno-unknown-pragmas -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc ConformalUtils.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
  icc -std=gnu++0x -O3 -openmp -o main main.cc Track.cc Hit.cc Matrix.cc KalmanUtils.cc Propagation.cc Simulation.cc buildtest.cc fittest.cc -I. `root-config --libs --cflags`
*/

#include <sstream>
#include <chrono>
#include <list>
#include <memory>

#include "Matrix.h"
#include "Event.h"
#include "Validation.h"
#include "BinInfoUtils.h"

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
  for (int l = 0; l < Config::nLayers; l++) {
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

  for (int l = 0; l < Config::nLayers; l++) {
    float r = (l+1)*Config::fRadialSpacing;
    float z = r / std::tan(2.0*std::atan(std::exp(-Config::fEtaDet))); // calculate z extent based on eta, r
    VUSolid* utub = new VUSolid(r, r+Config::fRadialExtent, -z, z, true, l + 1 == Config::nLayers);
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

// from mkFit
namespace
{
  std::string s_operation = "empty";
  std::string s_file_name = "simtracks.bin";
}

// also from mkfit
typedef std::list<std::string> lStr_t;
typedef lStr_t::iterator       lStr_i;

void next_arg_or_die(lStr_t& args, lStr_i& i, bool allow_single_minus=false)
{
  lStr_i j = i;
  if (++j == args.end() ||
      ((*j)[0] == '-' && ! (*j == "-" && allow_single_minus)))
  {
    std::cerr <<"Error: option "<< *i <<" requires an argument.\n";
    exit(1);
  }
  i = j;
}

int main(int argc, const char* argv[])
{

#ifdef TBB
  auto nThread(tbb::task_scheduler_init::default_num_threads());
#else
  auto nThread = 1;
#endif

  // following mkFit on argv
  lStr_t mArgs; 
  for (int i = 1; i < argc; ++i)
  {
    mArgs.push_back(argv[i]);
  }

  lStr_i i  = mArgs.begin();
  while (i != mArgs.end())
  {
    lStr_i start = i;

    if (*i == "-h" || *i == "-help" || *i == "--help")
    {
      printf(
        "Usage: %s [options]\n"
        "Options:\n"
        "  --num-events    <num>    number of events to run over (def: %d)\n"
        "  --num-tracks    <num>    number of tracks to generate for each event (def: %d)\n"
	"  --num-thr       <num>    number of threads used for TBB  (def: %d)\n"
	"  --root-val               bool to enable normal validation (eff, FR, DR) (def: %s)\n"
	"  --inc-shorts             include short reco tracks into FR (def: %s)\n"
	"  --cf-seeding             bool to enable CF in MC seeding (def: %s)\n"
	"  --read                   read input simtracks file (def: false)\n"
	"  --file-name              file name for write/read (def: %s)\n"
	"  --cmssw-seeds            take seeds from CMSSW (def: %s)\n"
        ,
        argv[0],
        Config::nEvents,
        Config::nTracks,
        nThread, 
	(Config::root_val ? "true" : "false"),
      	(Config::inclusiveShorts ? "true" : "false"),
	(Config::cf_seeding ? "true" : "false"),
	s_file_name.c_str(),
	(Config::seedInput == cmsswSeeds ? "true" : "false")
      );
      exit(0);
    }
    else if (*i == "--num-events")
    {
      next_arg_or_die(mArgs, i);
      Config::nEvents = atoi(i->c_str());
    }
    else if (*i == "--num-tracks")
    {
      next_arg_or_die(mArgs, i);
      Config::nTracks = atoi(i->c_str());
    }
    else if (*i == "--num-thr")
    {
      next_arg_or_die(mArgs, i);
      nThread = atoi(i->c_str());
    }
    else if (*i == "--root-val")
    {
      Config::root_val = true; Config::fit_val = false;
    }
    else if (*i == "--inc-shorts")
    {
      Config::inclusiveShorts = true;
    }
    else if (*i == "--cf-seeding")
    {
      Config::cf_seeding = true;
    }
    else if (*i == "--read")
    {
      s_operation = "read";
    }
    else if (*i == "--file-name")
    {
      next_arg_or_die(mArgs, i);
      s_file_name = *i;
    }
    else if(*i == "--cmssw-seeds")
    {
      Config::seedInput = cmsswSeeds;
    }
    else
    {
      fprintf(stderr, "Error: Unknown option/argument '%s'.\n", i->c_str());
      exit(1);
    }
    mArgs.erase(start, ++i);
  }

  Geometry geom;
  initGeom(geom);
  std::unique_ptr<Validation> val(Validation::make_validation("valtree.root"));

  for ( int i = 0; i < Config::nLayers; ++i ) {
    std::cout << "Layer = " << i << ", Radius = " << geom.Radius(i) << std::endl;
  }

  std::vector<tick> ticks(6);
  std::vector<unsigned int> tracks(4);
#ifdef TBB
  std::cout << "Initializing with " << nThread << " threads." << std::endl;
  tbb::task_scheduler_init tasks(nThread);
#endif

  DataFile data_file;
  if (s_operation == "read")
  {
    Config::nEvents = data_file.OpenRead(s_file_name);
  }

  for (int evt=0; evt<Config::nEvents; ++evt) {
    Event ev(geom, *val, evt, nThread);
    std::cout << "EVENT #"<< ev.evtID() << std::endl;

    timepoint t0(now());
    if (s_operation != "read")
    {
      ev.Simulate();
    }
    else {
      ev.read_in(data_file);
    }
    
    // phi-eta partitioning map: vector of vector of vectors of std::pairs. 
    // vec[nLayers][nEtaBins][nPhiBins]
    BinInfoMap segmentMap;
    
    /*simulate time*/        ticks[0] += delta(t0);
    ev.Segment(segmentMap);  ticks[1] += delta(t0);
    ev.Seed(segmentMap);     ticks[2] += delta(t0);
    ev.Find(segmentMap);     ticks[3] += delta(t0);
    ev.Fit();                ticks[4] += delta(t0);
    ev.Validate();           ticks[5] += delta(t0);

    std::cout << "sim: " << ev.simTracks_.size() << " seed: " << ev.seedTracks_.size() << " found: " 
	      << ev.candidateTracks_.size() << " fit: " << ev.fitTracks_.size() << std::endl;
    tracks[0] += ev.simTracks_.size();
    tracks[1] += ev.seedTracks_.size();
    tracks[2] += ev.candidateTracks_.size();
    tracks[3] += ev.fitTracks_.size();
    std::cout << "Built tracks" << std::endl;
    ev.PrintStats(ev.candidateTracks_, ev.candidateTracksExtra_);
    std::cout << "Fit tracks" << std::endl;
    ev.PrintStats(ev.fitTracks_, ev.fitTracksExtra_);
  }

  if (s_operation == "read")
  {
    data_file.Close();
  }

  std::vector<double> time(ticks.size());
  for (unsigned int i = 0; i < ticks.size(); ++i){
    time[i]=ticks[i].count();
  }

  val->fillConfigTree();
  val->saveTTrees(); 

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
