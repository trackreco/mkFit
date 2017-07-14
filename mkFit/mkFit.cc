#include "Matriplex/MatriplexCommon.h"

#include "fittestMPlex.h"
#include "buildtestMPlex.h"

#include "MkBuilder.h"
#include "MkFitter.h"

#include "Config.h"

#include "Timing.h"

#include <limits>
#include <list>
#include <sstream>
#include <memory>

#include "Event.h"

#include "MaterialEffects.h"

#ifndef NO_ROOT
#include "Validation.h"
#endif

#ifdef USE_CUDA
#include "FitterCU.h"
#include "gpu_utils.h"
#endif

#include <cstdlib>
//#define DEBUG
#include "Debug.h"

#include <tbb/task_scheduler_init.h>

#if defined(USE_VTUNE_PAUSE)
#include "ittnotify.h"
#endif

//==============================================================================

void initGeom(Geometry& geom)
{
  std::cout << "Constructing SimpleGeometry Cylinder geometry" << std::endl;

  // NB: we currently assume that each node is a layer, and that layers
  // are added starting from the center
  // NB: z is just a dummy variable, VUSolid is actually infinite in size.  *** Therefore, set it to the eta of simulation ***

  TrackerInfo::ExecTrackerInfoCreatorPlugin(Config::geomPlugin, Config::TrkInfo);

#ifndef WITH_USOLIDS
    geom.BuildFromTrackerInfo(Config::TrkInfo);
#else
    fprintf(stderr, "TrackerInfo only supports SimpleGeometry, not USolids.\n");
    exit(1);
#endif

  /*
  if ( ! Config::endcapTest && ! Config::useCMSGeom)
  {
    // This is the new standalone case -- Cylindrical Cow with Lids
    //Create_TrackerInfo(Config::TrkInfo);
#ifndef WITH_USOLIDS
    geom.BuildFromTrackerInfo(Config::TrkInfo);
#else
    fprintf(stderr, "Cylindrical Cow with Lids only supported for SimpleGeometry.\n");
    exit(1);
#endif
  }
  else
  {
    float eta = 2.0; // can tune this to whatever geometry required (one can make this layer dependent as well)
    for (int l = 0; l < Config::nLayers; l++)
    {
      if (Config::endcapTest)
      {
        float z = Config::useCMSGeom ? Config::cmsAvgZs[l] : (l+1)*10.;//Config::fLongitSpacing
        float rmin = Config::useCMSGeom ? Config::cmsDiskMinRs[l] : 0;
        float rmax = Config::useCMSGeom ? Config::cmsDiskMaxRs[l] : 0;
        // XXXX MT: Do we need endcap layer "thickness" for cmssw at all? Use equiv of fRadialExtent.
        // Do we even need geometry for cmssw?
        float dz = 0.005;
        VUSolid* utub = new VUSolid(rmin, rmax, z - dz, z + dz, false, l + 1 == Config::nLayers);
        geom.AddLayer(utub, rmin, z);
      }
      else
      {
        float r = Config::useCMSGeom ? Config::cmsAvgRads[l] : (l+1)*Config::fRadialSpacing;
        float z = r / std::tan(2.0*std::atan(std::exp(-eta))); // calculate z extent based on eta, r
        VUSolid* utub = new VUSolid(r, r+Config::fRadialExtent, -z, z, true, l + 1 == Config::nLayers);
        geom.AddLayer(utub, r, z);
      }
    }
  }
  */
}

namespace
{
  int   g_start_event   = 1;
  bool  g_run_fit_std   = false;

  bool  g_run_build_all = true;
  bool  g_run_build_bh  = false;
  bool  g_run_build_std = false;
  bool  g_run_build_ce  = false;

  std::string g_operation = "simulate_and_process";;
  std::string g_file_name = "simtracks.bin";
  std::string g_input_file = "";

  const char* b2a(bool b) { return b ? "true" : "false"; }
}

//==============================================================================

void read_and_save_tracks()
{
  DataFile in;
  const int Nevents = in.OpenRead(g_input_file, true);

  DataFile out;
  out.OpenWrite(g_file_name, Nevents);

  printf("writing %i events\n", Nevents);

  Event ev(0);
  for (int evt = 0; evt < Nevents; ++evt)
  {
    ev.Reset(evt);
    ev.read_in(in);
    ev.write_out(out);
  }

  out.Close();
  in .Close();
}

//==============================================================================

void generate_and_save_tracks()
{
  const int Nevents = Config::nEvents;

  Geometry geom;
  initGeom(geom);
  std::unique_ptr<Validation> val(Validation::make_validation("empty.root"));

  int extra_sections = 0;
  if (Config::root_val || Config::fit_val)
  {
    extra_sections |= DataFile::ES_SimTrackStates;
  }

  DataFile data_file;
  data_file.OpenWrite(g_file_name, Nevents, extra_sections);

  printf("writing %i events\n", Nevents);

  tbb::task_scheduler_init tbb_init(Config::numThreadsSimulation);

  Event ev(geom, *val, 0);
  for (int evt = 0; evt < Nevents; ++evt)
  {
    ev.Reset(evt);
    ev.Simulate();

#ifdef DEBUG
    for (int itrack = 0; itrack < ev.simTracks_.size(); itrack++)
    {
      const auto& track = ev.simTracks_[itrack];
      int mcTrackId = track.label();
      dprint("track: " << mcTrackId << " (" << itrack << ")");
      for (int ihit = 0; ihit < track.nTotalHits(); ihit++)
      {
	int idx = track.getHitIdx(ihit); int lyr = track.getHitLyr(ihit);
	int mcHitID = ev.layerHits_[lyr][idx].mcHitID(); int mcTrackID = ev.simHitsInfo_[mcHitID].mcTrackID();
	float tsr = ev.simTrackStates_[mcHitID].posR();	float hitr = ev.layerHits_[lyr][idx].r();
	float tsz = ev.simTrackStates_[mcHitID].z();    float hitz = ev.layerHits_[lyr][idx].z();
	dprint("       " << mcTrackID << " (mcHitID: " << mcHitID << " ihit: " << ihit << " idx: " << idx << " lyr: "
	       << lyr << " tsr: " << tsr << " hitr: " << hitr << " tsz: " << tsz << " hitz: " << hitz << ")");
      }
    }
#endif

    ev.write_out(data_file);
  }

  data_file.Close();
}

//==============================================================================

void test_standard()
{
  printf("Running test_standard(), operation=\"%s\"\n", g_operation.c_str());
  printf("  vusize=%d, num_th_sim=%d, num_th_finder=%d\n",
         MPT_SIZE, Config::numThreadsSimulation, Config::numThreadsFinder);
  printf("  sizeof(Track)=%zu, sizeof(Hit)=%zu, sizeof(SVector3)=%zu, sizeof(SMatrixSym33)=%zu, sizeof(MCHitInfo)=%zu\n",
         sizeof(Track), sizeof(Hit), sizeof(SVector3), sizeof(SMatrixSym33), sizeof(MCHitInfo));

  if (Config::useCMSGeom)     printf ("- using CMS-like geometry\n");
  if (Config::readCmsswSeeds) printf ("- reading seeds from file\n");
  if (Config::endcapTest)     printf ("- endcap test enabled (to be nixed)\n");

  if (g_operation == "write") {
    generate_and_save_tracks();
    return;
  }

  if (g_operation == "convert") {
    read_and_save_tracks();
    return;
  }

  Geometry geom;
  initGeom(geom);

  DataFile data_file;
  if (g_operation == "read")
  {
    int evs_in_file   = data_file.OpenRead(g_file_name);
    int evs_available = evs_in_file - g_start_event + 1;
    if (Config::nEvents == -1)
    {
      Config::nEvents = evs_available;
    }
    else if (Config::nEvents > evs_available)
    {
      printf("Requested number of events %d, only %d available.\n",
             Config::nEvents, evs_available);
      Config::nEvents = evs_available;
    }

    if (g_start_event > 1)
    {
      data_file.SkipNEvents(g_start_event - 1);
    }
  }

  if (Config::useCMSGeom) fillZRgridME();

  const int NT = 4;
  double t_sum[NT] = {0};
  double t_skip[NT] = {0};
  double time = dtime();

#if USE_CUDA
  tbb::task_scheduler_init tbb_init(Config::numThreadsFinder);
  //tbb::task_scheduler_init tbb_init(tbb::task_scheduler_init::automatic);

  //omp_set_num_threads(Config::numThreadsFinder);
  // fittest time. Sum of all events. In case of multiple events
  // being run simultaneously in different streams this time will
  // be larger than the elapsed time.

  std::vector<Event> events;
  std::vector<Validation> validations(Config::nEvents);

  events.reserve(Config::nEvents);
  // simulations are all performed before the fitting loop.
  // It is mandatory in order to see the benefits of running
  // multiple streams.
  for (int evt = 1; evt <= Config::nEvents; ++evt) {
    printf("Simulating event %d\n", evt);
    events.emplace_back(geom, val, evt);
    events.back().Simulate();
    dprint("Event #" << events.back().evtID() << " simtracks " << events.back().simTracks_.size() << " layerhits " << events.back().layerHits_.size());
  }

  // The first call to a GPU function always take a very long time.
  // Everything needs to be initialized.
  // Nothing is done in this function, except calling an harmless
  // CUDA function. These function can be changed to another one
  // if it becomes important to time (e.g. if you want the profiler to
  // tell you how much time is spend running cudaDeviceSynchronize(),
  // use another function).
  separate_first_call_for_meaningful_profiling_numbers();

  if (g_run_fit_std) runAllEventsFittingTestPlexGPU(events);

  if (g_run_build_all || g_run_build_bh) {
    double total_best_hit_time = 0.;
    total_best_hit_time = runAllBuildingTestPlexBestHitGPU(events);
    std::cout << "Total best hit time (GPU): " << total_best_hit_time << std::endl;
  }
#else
  std::atomic<int> nevt{g_start_event};
  std::atomic<int> seedstot{0}, simtrackstot{0};

  std::vector<std::unique_ptr<Event>>      evs(Config::numThreadsEvents);
  std::vector<std::unique_ptr<Validation>> vals(Config::numThreadsEvents);
  std::vector<std::unique_ptr<MkBuilder>>  mkbs(Config::numThreadsEvents);
  std::vector<std::shared_ptr<FILE>>       fps;
  fps.reserve(Config::numThreadsEvents);

  const std::string valfile("valtree");

  for (int i = 0; i < Config::numThreadsEvents; ++i) {
    std::ostringstream serial;
    if (Config::numThreadsEvents > 1) { serial << "_" << i; }
    vals[i].reset(Validation::make_validation(valfile + serial.str() + ".root"));
    mkbs[i].reset(MkBuilder::make_builder());
    evs[i].reset(new Event(geom, *vals[i], 0));
    if (g_operation == "read") {
      fps.emplace_back(fopen(g_file_name.c_str(), "r"), [](FILE* fp) { if (fp) fclose(fp); });
    }
  }

  tbb::task_scheduler_init tbb_init(Config::numThreadsFinder);

  dprint("parallel_for step size " << (Config::nEvents+Config::numThreadsEvents-1)/Config::numThreadsEvents);

  time = dtime();

  int events_per_thread = (Config::nEvents+Config::numThreadsEvents-1)/Config::numThreadsEvents;
  tbb::parallel_for(tbb::blocked_range<int>(0, Config::numThreadsEvents, 1),
    [&](const tbb::blocked_range<int>& threads)
  {
    int thisthread = threads.begin();

    assert(threads.begin() == threads.end()-1 && thisthread < Config::numThreadsEvents);

    std::vector<Track> plex_tracks;
    auto& ev     = *evs[thisthread].get();
    auto& mkb    = *mkbs[thisthread].get();
    auto  fp     =  fps[thisthread].get();

    int evstart = thisthread*events_per_thread;
    int evend   = std::min(Config::nEvents, evstart+events_per_thread);

    dprint("thisthread " << thisthread << " events " << Config::nEvents << " events/thread " << events_per_thread
                         << " range " << evstart << ":" << evend);

    for (int evt = evstart; evt < evend; ++evt)
    {
      ev.Reset(nevt++);

      if (!Config::silent)
      {
        std::lock_guard<std::mutex> printlock(Event::printmutex);
        printf("\n");
        printf("Processing event %d\n", ev.evtID());
      }

      if (g_operation == "read")
      {
        ev.read_in(data_file, fp);
      }
      else
      {
        ev.Simulate();
      }

      plex_tracks.resize(ev.simTracks_.size());

      double t_best[NT] = {0}, t_cur[NT];
      simtrackstot += ev.simTracks_.size();
      seedstot     += ev.seedTracks_.size();

      for (int b = 0; b < Config::finderReportBestOutOfN; ++b)
      {
  #ifndef USE_CUDA
        t_cur[0] = (g_run_fit_std) ? runFittingTestPlex(ev, plex_tracks) : 0;
  #else
        FitterCU<float> cuFitter(NN);
        cuFitter.allocateDevice();
        t_cur[0] = (g_run_fit_std) ? runFittingTestPlexGPU(cuFitter, ev, plex_tracks) : 0;
        cuFitter.freeDevice();
  #endif
        t_cur[1] = (g_run_build_all || g_run_build_bh)  ? runBuildingTestPlexBestHit(ev, mkb) : 0;
        t_cur[2] = (g_run_build_all || g_run_build_std) ? runBuildingTestPlexStandard(ev, mkb) : 0;
        t_cur[3] = (g_run_build_all || g_run_build_ce)  ? runBuildingTestPlexCloneEngine(ev, mkb) : 0;

        for (int i = 0; i < NT; ++i) t_best[i] = (b == 0) ? t_cur[i] : std::min(t_cur[i], t_best[i]);

        if (!Config::silent) {
          std::lock_guard<std::mutex> printlock(Event::printmutex);
          if (Config::finderReportBestOutOfN > 1)
          {
            printf("----------------------------------------------------------------\n");
            printf("Best-of-times:");
            for (int i = 0; i < NT; ++i) printf("  %.5f/%.5f", t_cur[i], t_best[i]);
            printf("\n");
          }
          printf("----------------------------------------------------------------\n");
        }
      }

      if (!Config::silent) {
        std::lock_guard<std::mutex> printlock(Event::printmutex);
        printf("Matriplex fit = %.5f  --- Build  BHMX = %.5f  STDMX = %.5f  CEMX = %.5f\n",
               t_best[0], t_best[1], t_best[2], t_best[3]);
      }

      // not protected by a mutex, may be inacccurate for multiple events in flight;
      // probably should convert to a scaled long so can use std::atomic<Integral>
      for (int i = 0; i < NT; ++i) t_sum[i] += t_best[i];
      if (evt > 0) for (int i = 0; i < NT; ++i) t_skip[i] += t_best[i];
    }
  }, tbb::simple_partitioner());

#endif
  time = dtime() - time;

  printf("\n");
  printf("================================================================\n");
  printf("=== TOTAL for %d events\n", Config::nEvents);
  printf("================================================================\n");

  printf("Total Matriplex fit = %.5f  --- Build  BHMX = %.5f  STDMX = %.5f  CEMX = %.5f\n",
         t_sum[0], t_sum[1], t_sum[2], t_sum[3]);
  printf("Total event > 1 fit = %.5f  --- Build  BHMX = %.5f  STDMX = %.5f  CEMX = %.5f\n",
         t_skip[0], t_skip[1], t_skip[2], t_skip[3]);
  printf("Total event loop time %.5f simtracks %d seedtracks %d\n", time, simtrackstot.load(), seedstot.load());
  //fflush(stdout);

  if (g_operation == "read")
  {
    data_file.Close();
  }

  for (auto& val : vals) {
    val->fillConfigTree();
    val->saveTTrees();
  }
}

//==============================================================================
// Command line argument parsing
//==============================================================================

typedef std::list<std::string> lStr_t;
typedef lStr_t::iterator       lStr_i;

bool has_suffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void next_arg_or_die(lStr_t& args, lStr_i& i, bool allow_single_minus=false)
{
  lStr_i j = i;
  if (++j == args.end() || has_suffix(*j, ".C") ||
      ((*j)[0] == '-' && ! (*j == "-" && allow_single_minus)))
  {
    std::cerr <<"Error: option "<< *i <<" requires an argument.\n";
    exit(1);
  }
  i = j;
}

//==============================================================================
// main
//==============================================================================

int main(int argc, const char *argv[])
{
#ifdef USE_VTUNE_PAUSE
  __itt_pause();
#endif

  assert (sizeof(Track::Status) == 4 && "To make sure this is true for icc and gcc<6 when mixing bools/ints in bitfields.");

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
        "  --geom          <str>    geometry plugin to use (def: %s)\n"
        "  --num-events    <num>    number of events to run over (def: %d)\n"
        "  --num-tracks    <num>    number of tracks to generate for each event (def: %d)\n"
        "  --num-thr-sim   <num>    number of threads for simulation (def: %d)\n"
        "  --num-thr       <num>    number of threads for track finding (def: %d)\n"
        "  --num-thr-ev    <num>    number of threads to run the event loop\n"
        "  --fit-std                run standard fitting test (def: false)\n"
        "  --fit-std-only           run only standard fitting test (def: false)\n"
        "  --chi2cut       <num>    chi2 cut used in building test (def: %.1f)\n"
        "  --build-bh               run best-hit building test (def: false)\n"
        "  --build-std              run standard combinatorial building test (def: false)\n"
        "  --build-ce               run clone engine combinatorial building test (def: false)\n"
        "  --seeds-per-task         number of seeds to process in a tbb task (def: %d)\n"
        "  --best-out-of   <num>    run track finding num times, report best time (def: %d)\n"
        "  --ext-rec-tracks         read external rec trakcs if available (def: %s)\n"
        "  --cmssw-seeds            take seeds from CMSSW (def: %s)\n"
        "  --find-seeds             run road search seeding [CF enabled by default] (def: %s)\n"
        "  --hits-per-task <num>    number of layer1 hits per task in finding seeds (def: %i)\n"
        "  --endcap-test            test endcap tracking (def: %s)\n"
        "  --cf-seeding             enable CF in seeding (def: %s)\n"
        "  --cf-fitting             enable CF in fitting (def: %s)\n"
        "  --root-val               enable ROOT based validation for building [eff, FR, DR] (def: %s)\n"
      	"  --fit-val                enable ROOT based validation for fitting (def: %s)\n"
        "  --inc-shorts             include short reco tracks into FR (def: %s)\n"
        "  --silent                 suppress printouts inside event loop (def: %s)\n"
        "  --write                  write simulation to file and exit\n"
        "  --read                   read simulation from file\n"
        "  --start-event   <num>    event number to start at when reading from a file (def: %d)\n"
        "  --file-name              file name for write/read (def: %s)\n"
        "  --input-file             file name for reading when converting formats (def: %s)\n"
        "GPU specific options: \n"
        "  --num-thr-reorg <num>    number of threads to run the hits reorganization\n"
        ,
        argv[0],
        Config::geomPlugin.c_str(),
        Config::nEvents,
        Config::nTracks,
        Config::numThreadsSimulation, Config::numThreadsFinder,
        Config::chi2Cut,
        Config::numSeedsPerTask,
        Config::finderReportBestOutOfN,
        b2a(Config::readCmsswSeeds),
        b2a(Config::findSeeds),
      	Config::numHitsPerTask,
        b2a(Config::endcapTest),
        b2a(Config::cf_seeding),
        b2a(Config::cf_fitting),
        b2a(Config::root_val),
        b2a(Config::fit_val),
	b2a(Config::inclusiveShorts),
        b2a(Config::silent),
        g_start_event,
      	g_file_name.c_str(),
      	g_input_file.c_str()
      );
      exit(0);
    }
    else if (*i == "--geom")
    {
      next_arg_or_die(mArgs, i);
      Config::geomPlugin = *i;
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
    else if (*i == "--num-thr-sim")
    {
      next_arg_or_die(mArgs, i);
      Config::numThreadsSimulation = atoi(i->c_str());
    }
    else if (*i == "--num-thr")
    {
      next_arg_or_die(mArgs, i);
      Config::numThreadsFinder = atoi(i->c_str());
    }
    else if(*i == "--fit-std")
    {
      g_run_fit_std = true;
    }
    else if(*i == "--fit-std-only")
    {
      g_run_fit_std = true;
      g_run_build_all = false; g_run_build_bh = false; g_run_build_std = false; g_run_build_ce = false;
    }
    else if (*i == "--chi2cut")
    {
      next_arg_or_die(mArgs, i);
      Config::chi2Cut = atof(i->c_str());
    }
    else if(*i == "--build-bh")
    {
      g_run_build_all = false; g_run_build_bh = true; g_run_build_std = false; g_run_build_ce = false;
    }
    else if(*i == "--build-std")
    {
      g_run_build_all = false; g_run_build_bh = false; g_run_build_std = true; g_run_build_ce = false;
    }
    else if(*i == "--build-ce")
    {
      g_run_build_all = false; g_run_build_bh = false; g_run_build_std = false; g_run_build_ce = true;
    }
    else if (*i == "--seeds-per-task")
    {
      next_arg_or_die(mArgs, i);
      Config::numSeedsPerTask = atoi(i->c_str());
    }
    else if(*i == "--best-out-of")
    {
      next_arg_or_die(mArgs, i);
      Config::finderReportBestOutOfN = atoi(i->c_str());
    }
    else if(*i == "--ext-rec-tracks")
    {
      Config::readExtRecTracks = true;
    }
    else if(*i == "--cmssw-seeds")
    {
      Config::readCmsswSeeds = true;
    }
    else if(*i == "--find-seeds")
    {
      Config::findSeeds = true; Config::cf_seeding = true;
    }
    else if (*i == "--hits-per-task")
    {
      next_arg_or_die(mArgs, i);
      Config::numHitsPerTask = atoi(i->c_str());
    }
    else if(*i == "--endcap-test")
    {
      Config::endcapTest = true; Config::nlayers_per_seed = 2; // default is 3 for barrel
    }
    else if (*i == "--cf-seeding")
    {
      Config::cf_seeding = true;
    }
    else if (*i == "--cf-fitting")
    {
      Config::cf_fitting = true;
    }
    else if (*i == "--root-val")
    {
      Config::root_val = true; Config::fit_val = false;
    }
    else if (*i == "--fit-val")
    {
      Config::root_val = false; Config::fit_val = true;
    }
    else if (*i == "--inc-shorts")
    {
      Config::inclusiveShorts = true;
    }
    else if (*i == "--num-thr-ev")
    {
      next_arg_or_die(mArgs, i);
      Config::numThreadsEvents = atoi(i->c_str());
    }
    else if (*i == "--num-thr-reorg")
    {
      next_arg_or_die(mArgs, i);
      Config::numThreadsReorg = atoi(i->c_str());
    }
    else if(*i == "--write")
    {
      g_operation = "write";
    }
    else if(*i == "--read")
    {
      g_operation = "read";
      Config::nEvents = -1;
    }
    else if (*i == "--start-event")
    {
      next_arg_or_die(mArgs, i);
      g_start_event = atoi(i->c_str());
    }
    else if(*i == "--file-name")
    {
      next_arg_or_die(mArgs, i);
      g_file_name = *i;
    }
    else if(*i == "--input-file")
    {
      next_arg_or_die(mArgs, i);
      g_operation = "convert";
      g_input_file = *i;
    }
    else if(*i == "--silent")
    {
      Config::silent = true;
    }
    else
    {
      fprintf(stderr, "Error: Unknown option/argument '%s'.\n", i->c_str());
      exit(1);
    }

    mArgs.erase(start, ++i);
  }

  Config::RecalculateDependentConstants();

  printf ("Running with n_threads=%d, best_out_of=%d\n",
          Config::numThreadsFinder, Config::finderReportBestOutOfN);

  test_standard();

  return 0;
}
