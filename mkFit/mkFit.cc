#include "Matriplex/MatriplexCommon.h"

#include "fittestMPlex.h"
#include "buildtestMPlex.h"

#include "MkFitter.h"

#include "Config.h"

#include "Timing.h"

#include <limits>

#include "Event.h"

#ifndef NO_ROOT
#include "TTreeValidation.h"
#endif

#ifdef USE_CUDA
#include "FitterCU.h"
#endif

#include <omp.h>

#include <tbb/task_scheduler_init.h>

#if defined(USE_VTUNE_PAUSE)
#include "ittnotify.h"
#endif

//==============================================================================

std::vector<Track> plex_tracks;

void initGeom(Geometry& geom)
{
  std::cout << "Constructing SimpleGeometry Cylinder geometry" << std::endl;

  // NB: we currently assume that each node is a layer, and that layers
  // are added starting from the center
  // NB: z is just a dummy variable, VUSolid is actually infinite in size.  *** Therefore, set it to the eta of simulation ***
  float eta = 2.0; // can tune this to whatever geometry required (one can make this layer dependent as well)
  for (int l = 0; l < Config::nLayers; l++) {
    float r = Config::useCMSGeom ? Config::cmsAvgRads[l] : (l+1)*Config::fRadialSpacing;
    VUSolid* utub = new VUSolid(r, r+Config::fRadialExtent);
    float z = r / std::tan(2.0*std::atan(std::exp(-eta))); // calculate z extent based on eta, r
    geom.AddLayer(utub, r, z);
  }
}

namespace
{
  FILE *g_file = 0;
  int   g_file_num_ev = 0;
  int   g_file_cur_ev = 0;

  bool  g_run_fit_std   = false;

  bool  g_run_build_all = true;
  bool  g_run_build_bh  = false;
  bool  g_run_build_std = false;
  bool  g_run_build_ce  = false;
  bool  g_run_build_tbb = false;

  std::string g_operation = "simulate_and_process";;
  std::string g_file_name = "simtracks.bin";
}

void generate_and_save_tracks()
{
  FILE *fp = fopen(g_file_name.c_str(), "w");

  int Ntracks = Config::nTracks;

  int Nevents = Config::nEvents;

  Geometry geom;
  initGeom(geom);
  Validation val;

  fwrite(&Nevents, sizeof(int), 1, fp);

  printf("writing %i events\n",Nevents);

  for (int evt = 0; evt < Nevents; ++evt)
  {
    Event ev(geom, val, evt);

    omp_set_num_threads(Config::numThreadsSimulation);

    ev.Simulate();
    ev.resetLayerHitMap(true);

    ev.write_out(fp);
  }

  fclose(fp);
}


int open_simtrack_file()
{
  g_file = fopen(g_file_name.c_str(), "r");

  assert (g_file != 0);

  fread(&g_file_num_ev, sizeof(int), 1, g_file);
  g_file_cur_ev = 0;

  printf("\nReading simulated tracks from \"%s\", %d events on file.\n\n",
         g_file_name.c_str(), g_file_num_ev);

  return g_file_num_ev;
}

void close_simtrack_file()
{
  fclose(g_file);
  g_file = 0;
  g_file_num_ev = 0;
  g_file_cur_ev = 0;
}

void test_standard()
{
  // ---- MT test eta bins
  // int nb, b1, b2;
  // for (float eta = -1.2; eta <= 1.2; eta += 0.01)
  // {
  //   nb = getBothEtaBins(eta, b1, b2);
  //   printf("eta=%6.2f  bin=%3d  bin1=%3d  bin2=%3d nb=%d\n",
  //          eta, getEtaBin(eta), b1, b2, nb);
  // }

  // return;
  // ---- end MT test

  printf("Running test_standard(), operation=\"%s\"\n", g_operation.c_str());
  printf("  vusize=%d, num_th_sim=%d, num_th_finder=%d\n",
         MPT_SIZE, Config::numThreadsSimulation, Config::numThreadsFinder);
  printf("  sizeof(Track)=%zu, sizeof(Hit)=%zu, sizeof(SVector3)=%zu, sizeof(SMatrixSym33)=%zu, sizeof(MCHitInfo)=%zu\n",
         sizeof(Track), sizeof(Hit), sizeof(SVector3), sizeof(SMatrixSym33), sizeof(MCHitInfo));
  if (Config::useCMSGeom) {
    printf ("Using CMS-like geometry ");
    if (Config::readCmsswSeeds) printf ("with CMSSW seeds \n");
    else printf ("with MC-truth seeds \n");
  } else printf ("Using 4-cm spacing geometry \n");

  if (g_operation == "write") {
    generate_and_save_tracks();
    return;
  }

  if (g_operation == "read")
  {
    Config::nEvents = open_simtrack_file();
    //Config::nEvents = 10;
  }

  Geometry geom;
  initGeom(geom);
#ifdef NO_ROOT
  Validation val;
#else 
  TTreeValidation val("valtree.root");
#endif
  
  const int NT = 5;
  double t_sum[NT] = {0};
  double t_skip[NT] = {0};

  EventTmp ev_tmp;

#if USE_CUDA
  // fittest time. Sum of all events. In case of multiple events
  // being run simultaneously in different streams this time will
  // be larger than the elapsed time.
  double s_tmp = 0.0;

  std::vector<Event> events;
  std::vector<Validation> validations(Config::nEvents);

  // simulations are all performed before the fitting loop.
  // It is mandatory in order to see the benefits of running
  // multiple streams.
  for (int evt = 1; evt <= Config::nEvents; ++evt) {
    printf("Simulating event %d\n", evt);
    Event ev(geom, val, evt);
    ev.Simulate();
    ev.resetLayerHitMap(true);
    
    events.push_back(ev);
  }

  // The first call to a GPU function always take a very long time.
  // Everything needs to be initialized.
  // Nothing is done in this function, except calling an harmless
  // CUDA function. These function can be changed to another one
  // if it becomes important to time (e.g. if you want the profiler to
  // tell you how much time is spend running cudaDeviceSynchronize(),
  // use another function). 
  separate_first_call_for_meaningful_profiling_numbers();
#if 0
  In principle, the warmup loop should not be required.
  The separate_first_call_for_meaningful_profiling_numbers() function
  should be enough.
  // Warmup loop
  for (int i = 0; i < 1; ++i) {
    FitterCU<float> cuFitter(NN);
    cuFitter.allocateDevice();
    Event &ev = events[0];
    std::vector<Track> plex_tracks_ev;
    plex_tracks_ev.resize(ev.simTracks_.size());

    if (g_run_fit_std) runFittingTestPlexGPU(cuFitter, ev, plex_tracks_ev);
    cuFitter.freeDevice();
  }
#endif

  // Reorgnanization (copyIn) can eventually be multithreaded.
  omp_set_nested(1);
      
  omp_set_num_threads(Config::numThreadsEvents);
  double total_gpu_time = dtime();
#pragma omp parallel reduction(+:s_tmp)
  {
  int numThreadsEvents = omp_get_num_threads();
  int thr_idx = omp_get_thread_num();

  // FitterCU is declared here to share allocations and deallocations
  // between the multiple events processed by a single thread.
  FitterCU<float> cuFitter(NN);
  cuFitter.allocateDevice();

    for (int evt = thr_idx+1; evt <= Config::nEvents; evt+= numThreadsEvents) {
      int idx = thr_idx;
      printf("==============================================================\n");
      printf("Processing event %d with thread %d\n", evt, idx);
      Event &ev = events[evt-1];
      std::vector<Track> plex_tracks_ev;
      plex_tracks_ev.resize(ev.simTracks_.size());
      double tmp = 0, tmp2bh = 0, tmp2 = 0, tmp2ce = 0;

      if (g_run_fit_std) tmp = runFittingTestPlexGPU(cuFitter, ev, plex_tracks_ev);

      printf("Matriplex fit = %.5f  -------------------------------------", tmp);
      printf("\n");
      s_tmp    += tmp;
#if 1  // 0 for timing, 1 for validation
      // Validation crashes for multiple threads.
      // It is something in relation to ROOT. Not sure what. 
      if (omp_get_num_threads() <= 1) {
        if (g_run_fit_std) {
          std::string tree_name = "validation-plex-" + std::to_string(evt) + ".root";
          make_validation_tree(tree_name.c_str(), ev.simTracks_, plex_tracks_ev);
        }
      }
#endif
    }
    cuFitter.freeDevice();
  }
  std::cerr << "###### Total GPU time: " << dtime() - total_gpu_time << " ######\n";

#else
  // MT: task_scheduler_init::automatic doesn't really work (segv!) + we don't
  // know what to do for non-tbb cases.
  // tbb::task_scheduler_init tbb_init(Config::numThreadsFinder != 0 ?
  //                                   Config::numThreadsFinder :
  //                                   tbb::task_scheduler_init::automatic);
  tbb::task_scheduler_init tbb_init(Config::numThreadsFinder);
  omp_set_num_threads(Config::numThreadsFinder);

  for (int evt = 1; evt <= Config::nEvents; ++evt)
  {
    printf("\n");
    printf("Processing event %d\n", evt);

    Event ev(geom, val, evt);

    if (g_operation == "read")
    {
      ev.read_in(g_file);
      ev.resetLayerHitMap(false);//hitIdx's in the sim tracks are already ok 
    }
    else
    {
      //Simulate() parallelism is via TBB, but comment out for now due to cost of
      //task_scheduler_init
      //tbb::task_scheduler_init tbb_init(Config::numThreadsSimulation);

      ev.Simulate();
      ev.resetLayerHitMap(true);
    }

    // if (evt!=2985) continue;

    plex_tracks.resize(ev.simTracks_.size());

    double t_best[NT] = {0}, t_cur[NT];

    for (int b = 0; b < Config::finderReportBestOutOfN; ++b)
    {
      t_cur[0] = (g_run_fit_std) ? runFittingTestPlex(ev, plex_tracks) : 0;
      t_cur[1] = (g_run_build_all || g_run_build_bh)  ? runBuildingTestPlexBestHit(ev) : 0;
      t_cur[2] = (g_run_build_all || g_run_build_std) ? runBuildingTestPlex(ev, ev_tmp) : 0;
      t_cur[3] = (g_run_build_all || g_run_build_ce)  ? runBuildingTestPlexCloneEngine(ev, ev_tmp) : 0;
      t_cur[4] = (g_run_build_all || g_run_build_tbb) ? runBuildingTestPlexTbb(ev, ev_tmp) : 0;

      for (int i = 0; i < NT; ++i) t_best[i] = (b == 0) ? t_cur[i] : std::min(t_cur[i], t_best[i]);

      if (Config::finderReportBestOutOfN > 1)
      {
        printf("----------------------------------------------------------------\n");
        printf("Best-of-times:");
        for (int i = 0; i < NT; ++i) printf("  %.5f/%.5f", t_cur[i], t_best[i]);
        printf("\n");
      }
      printf("----------------------------------------------------------------\n");
    }

    printf("Matriplex fit = %.5f  --- Build  BHMX = %.5f  MX = %.5f  CEMX = %.5f  TBBMX = %.5f\n",
           t_best[0], t_best[1], t_best[2], t_best[3], t_best[4]);

    for (int i = 0; i < NT; ++i) t_sum[i] += t_cur[i];
    if (evt > 1) for (int i = 0; i < NT; ++i) t_skip[i] += t_cur[i];

#ifndef NO_ROOT
    make_validation_tree("validation-plex.root", ev.simTracks_, plex_tracks);
#endif
    
  }
#endif
  printf("\n");
  printf("================================================================\n");
  printf("=== TOTAL for %d events\n", Config::nEvents);
  printf("================================================================\n");

  printf("Total Matriplex fit = %.5f  --- Build  BHMX = %.5f  MX = %.5f  CEMX = %.5f  TBBMX = %.5f\n",
         t_sum[0], t_sum[1], t_sum[2], t_sum[3], t_sum[4]);
  printf("Total event > 1 fit = %.5f  --- Build  BHMX = %.5f  MX = %.5f  CEMX = %.5f  TBBMX = %.5f\n",
         t_skip[0], t_skip[1], t_skip[2], t_skip[3], t_skip[4]);

  if (g_operation == "read")
  {
    close_simtrack_file();
  }

}

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

int main(int argc, const char *argv[])
{
#ifdef USE_VTUNE_PAUSE
  __itt_pause();
#endif

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
        "  --num-thr-sim   <num>    number of threads for simulation (def: %d)\n"
        "  --num-thr       <num>    number of threads for track finding (def: %d)\n"
        "                           extra cloning thread is spawned for each of them\n"
        "  --fit-std                run standard fitting test (def: false)\n"
        "  --fit-std-only           run only standard fitting test (def: false)\n"
        "  --build-bh               run best-hit building test (def: run all building tests)\n"
        "  --build-std              run standard building test\n"
        "  --build-ce               run clone-engine building test\n"
        "  --build-tbb              run tbb building test\n"
        "  --cloner-single-thread   do not spawn extra cloning thread (def: %s)\n"
        "  --seeds-per-task         number of seeds to process in a tbb task (def: %d)\n"
        "  --best-out-of   <num>    run track finding num times, report best time (def: %d)\n"
	"  --cms-geom               use cms-like geometry (def: %i)\n"
	"  --cmssw-seeds            take seeds from CMSSW (def: %i)\n"
	"  --cf-seeding             enable CF in seeding (def: %s)\n"
	"  --cf-fitting             enable CF in fitting (def: %s)\n"
	"  --write                  write simulation to file and exit\n"
	"  --read                   read simulation from file\n"
	"  --file-name              file name for write/read (def: %s)\n"
        "GPU specific options: \n"
        "  --num-thr-ev    <num>    number of threads to run the event loop\n"
        "  --num-thr-reorg <num>    number of threads to run the hits reorganization\n"
        ,
        argv[0],
        Config::nEvents,
        Config::nTracks,
        Config::numThreadsSimulation, Config::numThreadsFinder,
        Config::clonerUseSingleThread ? "true" : "false",
        Config::numSeedsPerTask,
        Config::finderReportBestOutOfN,
	Config::useCMSGeom,
	Config::readCmsswSeeds,
	Config::cf_seeding ? "true" : "false",
	Config::cf_fitting ? "true" : "false",
	g_file_name.c_str()
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
    else if(*i == "--build-bh")
    {
      g_run_build_all = false; g_run_build_bh = true;
    }
    else if(*i == "--build-std")
    {
      g_run_build_all = false; g_run_build_std = true;
    }
    else if(*i == "--build-ce")
    {
      g_run_build_all = false; g_run_build_ce = true;
    }
    else if(*i == "--build-tbb")
    {
      g_run_build_all = false; g_run_build_tbb = true;
    }
    else if(*i == "--cloner-single-thread")
    {
      Config::clonerUseSingleThread = true;
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
    else if(*i == "--cms-geom")
    {
      Config::useCMSGeom = true;
    }
    else if(*i == "--cmssw-seeds")
    {
      Config::readCmsswSeeds = true;
    }
    else if (*i == "--cf-seeding")
    {
      Config::cf_seeding = true;
    }
    else if (*i == "--cf-fitting")
    {
      Config::cf_fitting = true;
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
    }
    else if(*i == "--file-name")
    {
      next_arg_or_die(mArgs, i);
      g_file_name = *i;
    }
    else
    {
      fprintf(stderr, "Error: Unknown option/argument '%s'.\n", i->c_str());
      exit(1);
    }

    mArgs.erase(start, ++i);
  }

  Config::RecalculateDependentConstants();

  printf ("Running with n_threads=%d, cloner_single_thread=%d, best_out_of=%d\n",
          Config::numThreadsFinder, Config::clonerUseSingleThread, Config::finderReportBestOutOfN);

  test_standard();

  return 0;
}
