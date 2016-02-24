#include "Matriplex/MatriplexCommon.h"

#include "fittestMPlex.h"
#include "buildtestMPlex.h"

#include "MkFitter.h"

#include "Config.h"

#include "Timing.h"

#include <limits>

#include "Event.h"

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
  for (int l = 0; l < 10; l++) {
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

  std::string g_operation = "simulate_and_process";;
  std::string g_file_name = "simtracks.bin";
}

// take out the part for reading and writing the event
/*
void generate_and_save_tracks()
{
  FILE *fp = fopen(g_file_name.c_str(), "w");

  int Ntracks = Config::nTracks;

  int Nevents = Config::nEvents;

  Geometry geom;
  initGeom(geom);
  Validation val;

  fwrite(&Nevents, sizeof(int), 1, fp);

  for (int evt = 0; evt < Nevents; ++evt)
  {
    Event ev(geom, val, evt);

    ev.Simulate();
    ev.resetLayerHitMap();

    fwrite(&Ntracks, sizeof(int), 1, fp);

    for (int i = 0; i < Ntracks; ++i)
    {
      // ev.simTracks_[i].write_out(fp);
    }
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

int read_simtrack_event(std::vector<Track> &simtracks)
{
  int nt;

  fread(&nt, sizeof(int), 1, g_file);

  std::vector<Track> new_tracks(nt);
  simtracks.swap(new_tracks);

  for (int i = 0; i < nt; ++i)
  {
    // simtracks[i].read_in(g_file);
  }

  ++g_file_cur_ev;

  return nt;
}

void close_simtrack_file()
{
  fclose(g_file);
  g_file = 0;
  g_file_num_ev = 0;
  g_file_cur_ev = 0;
}
*/

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
  if (Config::useCMSGeom) printf ("Using CMS-like geometry \n");
  else printf ("Using 4-cm spacing geometry \n");


  if (g_operation == "read")
  {
    // Nevents = open_simtrack_file();
  }

  Geometry geom;
  initGeom(geom);
  Validation val;

  const int NT = 4;
  double t_sum[NT] = {0};

  EventTmp ev_tmp;

  tbb::task_scheduler_init tbb_init(Config::numThreadsFinder != 0 ?
                                    Config::numThreadsFinder :
                                    tbb::task_scheduler_init::automatic);

  for (int evt = 1; evt <= Config::nEvents; ++evt)
  {
    printf("\n");
    printf("Processing event %d\n", evt);

    Event ev(geom, val, evt);

    if (g_operation == "read")
    {
      // Ntracks = read_simtrack_event(ev.simTracks_);
    }
    else
    {
      omp_set_num_threads(Config::numThreadsSimulation);

      ev.Simulate();
      ev.resetLayerHitMap(true);

      omp_set_num_threads(Config::numThreadsFinder);
    }

    plex_tracks.resize(ev.simTracks_.size());

    double t_best[NT] = {0}, t_cur[NT];

    for (int b = 0; b < Config::finderReportBestOutOfN; ++b)
    {
      t_cur[0] = (g_run_fit_std) ? runFittingTestPlex(ev, plex_tracks) : 0;

      t_cur[1] = (g_run_build_all || g_run_build_bh)  ? runBuildingTestPlexBestHit(ev) : 0;

      t_cur[2] = (g_run_build_all || g_run_build_std) ? runBuildingTestPlex(ev, ev_tmp) : 0;

      t_cur[3] = (g_run_build_all || g_run_build_ce)  ? runBuildingTestPlexCloneEngine(ev, ev_tmp) : 0;

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

    printf("Matriplex fit = %.5f  --- Build  BHMX = %.5f  MX = %.5f  CEMX = %.5f\n",
           t_best[0], t_best[1], t_best[2], t_best[3]);

    for (int i = 0; i < NT; ++i) t_sum[i] += t_best[i];
  }
  printf("\n");
  printf("================================================================\n");
  printf("=== TOTAL for %d events\n", Config::nEvents);
  printf("================================================================\n");

  printf("Total Matriplex fit = %.5f  --- Build  BHMX = %.5f  MX = %.5f  CEMX = %.5f\n",
         t_sum[0], t_sum[1], t_sum[2], t_sum[3]);

  if (g_operation == "read")
  {
    // close_simtrack_file();
  }

#ifndef NO_ROOT
  make_validation_tree("validation-plex.root", ev.simTracks_, plex_tracks);
#endif
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
        "  --num-thr-sim   <num>    number of threads for simulation (def: %d)\n"
        "  --num-thr       <num>    number of threads for track finding (def: %d)\n"
        "                           extra cloning thread is spawned for each of them\n"
        "  --fit-std                run standard fitting test (def: false)\n"
        "  --fit-std-only           run only standard fitting test (def: false)\n"
        "  --build-bh               run best-hit building test (def: run all building tests)\n"
        "  --build-std              run standard building test\n"
        "  --build-ce               run clone-engine building test\n"
        "  --cloner-single-thread   do not spawn extra cloning thread (def: %s)\n"
        "  --seeds-per-task         number of seeds to process in a tbb task (def: %d)\n"
        "  --best-out-of   <num>    run track finding num times, report best time (def: %d)\n"
	"  --cms-geom               use cms-like geometry (def: %i)\n"
        ,
        argv[0],
        Config::numThreadsSimulation, Config::numThreadsFinder,
        Config::clonerUseSingleThread ? "true" : "false",
        Config::finderReportBestOutOfN,
        Config::numSeedsPerTask,
	Config::useCMSGeom
      );
      exit(0);
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
    else
    {
      fprintf(stderr, "Error: Unknown option/argument '%s'.\n", i->c_str());
      exit(1);
    }

    mArgs.erase(start, ++i);
  }

  printf ("Running with n_threads=%d, cloner_single_thread=%d, best_out_of=%d\n",
          Config::numThreadsFinder, Config::clonerUseSingleThread, Config::finderReportBestOutOfN);

  /*
  if (argc >= 2)
  {
    g_operation = argv[1];

    if (g_operation != "write" && g_operation != "read")
    {
      usage_and_die(argv[0]);
    }

    if (argc == 3)
    {
      g_file_name = argv[2];
    }

    if (argc > 3)
    {
      usage_and_die(argv[0]);
    }
  }

  if (g_operation == "write")
  {
    //fixme generate_and_save_tracks();
  }
  else
  {
    test_standard();
  }
  */

  test_standard();

  return 0;
}
