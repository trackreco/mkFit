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
#include "BuilderCU.h"
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
  if ( ! Config::useCMSGeom)
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
        float r = Config::useCMSGeom ? Config::cmsAvgRads[l] : (l+1)*Config::fRadialSpacing;
        float z = r / std::tan(2.0*std::atan(std::exp(-eta))); // calculate z extent based on eta, r
        VUSolid* utub = new VUSolid(r, r+Config::fRadialExtent, -z, z, true, l + 1 == Config::nLayers);
        geom.AddLayer(utub, r, z);
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
  bool  g_run_build_fv  = false;
  bool  g_seed_based    = false;

  std::string g_operation = "simulate_and_process";;
  std::string g_input_file = "";
  std::string g_output_file = "";

  seedOptsMap g_seed_opts;
  void init_seed_opts()
  {
    g_seed_opts["sim"]   = simSeeds;
    g_seed_opts["cmssw"] = cmsswSeeds;
    g_seed_opts["find"]  = findSeeds;
  }

  cleanOptsMap g_clean_opts;
  void init_clean_opts()
  {
    g_clean_opts["none"]     = noCleaning;
    g_clean_opts["n2"]       = cleanSeedsN2;
    g_clean_opts["pure"]     = cleanSeedsPure;
    g_clean_opts["badlabel"] = cleanSeedsBadLabel;
  }

  matchOptsMap g_match_opts;
  void init_match_opts()
  {
    g_match_opts["trkparam"] = trkParamBased;
    g_match_opts["hits"]     = hitBased;
    g_match_opts["label"]    = labelBased;
  }

  const char* b2a(bool b) { return b ? "true" : "false"; }
}

//==============================================================================

// Getters and setters of enum configs (from command line using anon. namespace above)

template <typename T, typename U> 
std::string getOpt(const T & c_opt, const U & g_opt_map)
{
  static const std::string empty("");

  for (const auto & g_opt_pair : g_opt_map)
  {
    if (g_opt_pair.second == c_opt) return g_opt_pair.first;
  }
  std::cerr << "No match for option " << c_opt << std::endl;
  return empty;
}

template <typename T, typename U>
void setOpt(const std::string & cmd_ln_str, T & c_opt, const U & g_opt_map, const std::string & ex_txt)
{
  if (g_opt_map.count(cmd_ln_str)) c_opt = g_opt_map.at(cmd_ln_str);
  else 
  {
    std::cerr << cmd_ln_str << " is not a valid " << ex_txt << " option!! Exiting..." << std::endl;
    exit(1);
  }
}

//==============================================================================

void read_and_save_tracks()
{
  DataFile in;
  const int Nevents = in.OpenRead(g_input_file, true);

  DataFile out;
  out.OpenWrite(g_output_file, Nevents);

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
  data_file.OpenWrite(g_output_file, Nevents, extra_sections);

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

  if (Config::useCMSGeom)              printf ("- using CMS-like geometry\n");
  if (Config::seedInput == cmsswSeeds) printf ("- reading seeds from file\n");

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
    int evs_in_file   = data_file.OpenRead(g_input_file);
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

  const int NT = 5;
  double t_sum[NT] = {0};
  double t_skip[NT] = {0};
  double time = dtime();

#if USE_CUDA_OLD
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

  MkBuilder::populate(g_run_build_all || g_run_build_fv);

  std::vector<std::unique_ptr<Event>>      evs(Config::numThreadsEvents);
  std::vector<std::unique_ptr<Validation>> vals(Config::numThreadsEvents);
  std::vector<std::unique_ptr<MkBuilder>>  mkbs(Config::numThreadsEvents);
  std::vector<std::shared_ptr<FILE>>       fps;
  fps.reserve(Config::numThreadsEvents);

#if USE_CUDA
  separate_first_call_for_meaningful_profiling_numbers();

  std::vector<std::unique_ptr<FitterCU<float>>> cuFitters(Config::numThreadsEvents);
  std::vector<std::unique_ptr<BuilderCU>> cuBuilders(Config::numThreadsEvents);
#endif

  const std::string valfile("valtree");

  for (int i = 0; i < Config::numThreadsEvents; ++i) {
    std::ostringstream serial;
    if (Config::numThreadsEvents > 1) { serial << "_" << i; }
    vals[i].reset(Validation::make_validation(valfile + serial.str() + ".root"));
    mkbs[i].reset(MkBuilder::make_builder());
    evs[i].reset(new Event(geom, *vals[i], 0));
    if (g_operation == "read") {
      fps.emplace_back(fopen(g_input_file.c_str(), "r"), [](FILE* fp) { if (fp) fclose(fp); });
    }
#if USE_CUDA
    constexpr int gplex_width = 10000;
    cuFitters[i].reset(new FitterCU<float>(gplex_width));
    cuFitters[i].get()->allocateDevice();
    cuFitters[i].get()->allocate_extra_addBestHit();
    cuFitters[i].get()->allocate_extra_combinatorial();
    cuFitters[i].get()->createStream();
    cuFitters[i].get()->setNumberTracks(gplex_width);

    cuBuilders[i].reset(new BuilderCU(cuFitters[i].get()));
    cuBuilders[i].get()->allocateGeometry(geom);
#endif
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

#if USE_CUDA
    auto& cuFitter = *cuFitters[thisthread].get();
    auto& cuBuilder = *cuBuilders[thisthread].get();
#endif

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
        t_cur[1] = (g_run_build_all || g_run_build_bh)  ? runBuildingTestPlexBestHit(ev, mkb) : 0;
        t_cur[3] = (g_run_build_all || g_run_build_ce)  ? runBuildingTestPlexCloneEngine(ev, mkb) : 0;
        t_cur[4] = (g_run_build_all || g_run_build_fv)  ? runBuildingTestPlexFV(ev, mkb) : 0;
  #else
        t_cur[0] = (g_run_fit_std) ? runFittingTestPlexGPU(cuFitter, ev, plex_tracks) : 0;
        t_cur[1] = (g_run_build_all || g_run_build_bh)  ? runBuildingTestPlexBestHitGPU(ev, mkb, cuBuilder) : 0;
        // XXXX MT note for Matthieu: ev_tmp no longer exists ----------------------------------v
        t_cur[3] = (g_run_build_all || g_run_build_ce)  ? runBuildingTestPlexCloneEngineGPU(ev, ev_tmp, mkb, cuBuilder, g_seed_based) : 0;
  #endif
        t_cur[2] = (g_run_build_all || g_run_build_std) ? runBuildingTestPlexStandard(ev, mkb) : 0;

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
        printf("Matriplex fit = %.5f  --- Build  BHMX = %.5f  STDMX = %.5f  CEMX = %.5f  FVMX = %.5f\n",
               t_best[0], t_best[1], t_best[2], t_best[3], t_best[4]);
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

  printf("Total Matriplex fit = %.5f  --- Build  BHMX = %.5f  STDMX = %.5f  CEMX = %.5f  FVMX = %.5f\n",
         t_sum[0], t_sum[1], t_sum[2], t_sum[3], t_sum[4]);
  printf("Total event > 1 fit = %.5f  --- Build  BHMX = %.5f  STDMX = %.5f  CEMX = %.5f  FVMX = %.5f\n",
         t_skip[0], t_skip[1], t_skip[2], t_skip[3], t_skip[4]);
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
#if USE_CUDA
  for (int i = 0; i < Config::numThreadsEvents; ++i) {
    cuFitters[i].get()->freeDevice();
    cuFitters[i].get()->free_extra_addBestHit();
    cuFitters[i].get()->free_extra_combinatorial();
    cuFitters[i].get()->destroyStream();
  }
#endif
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

  // init enum maps
  init_seed_opts();
  init_clean_opts();
  init_match_opts();

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
        "  --num-thr-ev    <num>    number of threads to run the event loop (def: %d)\n"
        "  --fit-std                run standard fitting test (def: false)\n"
        "  --fit-std-only           run only standard fitting test (def: false)\n"
        "  --chi2cut       <num>    chi2 cut used in building test (def: %.1f)\n"
        "  --build-bh               run best-hit building test (def: false)\n"
        "  --build-std              run standard combinatorial building test (def: false)\n"
        "  --build-ce               run clone engine combinatorial building test (def: false)\n"
        "  --build-fv               run full vector combinatorial building test (def: false)\n"
	"  --use-phiq-arr           use phi-Q arrays in select hit indices (def: %s)\n"
        "  --seeds-per-task         number of seeds to process in a tbb task (def: %d)\n"
        "  --hits-per-task <num>    number of layer1 hits per task in finding seeds (def: %i)\n"
        "  --best-out-of   <num>    run track finding num times, report best time (def: %d)\n"
        "  --seed-input    <str>    which seed collecion used for building (def: %s)\n"
        "  --seed-cleaning <str>    which seed cleaning to apply if using cmssw seeds (def: %s)\n" 
        "  --cf-seeding             enable CF in seeding (def: %s)\n"
        "  --cf-fitting             enable CF in fitting (def: %s)\n"
        "  --read-cmssw-tracks      read external cmssw reco tracks if available (def: %s)\n"
	"  --read-simtrack-states   read in simTrackStates for pulls in validation (def: %s)\n"
	"  --quality-val            enable printout validation for MkBuilder (def: %s)\n"
        "  --root-val               enable ROOT based validation for building [eff, FR, DR] (def: %s)\n"
        "  --cmssw-val              enable special CMSSW ROOT based validation for building [eff] (def: %s)\n"
      	"  --fit-val                enable ROOT based validation for fitting (def: %s)\n"
        "  --inc-shorts             include short reco tracks into FR (def: %s)\n"
	"  --cmssw-matching <str>   which cmssw track matching routine to use if doing special cmssw validation, candidate tracks only (def: %s)\n" 
        "  --hit-match              apply hit matching criteria for track param matching in cmssw validation (def: %s)\n"
	"  --dump-for-plots         printouts for plots from logs (def: %s)\n"
        "  --silent                 suppress printouts inside event loop (def: %s)\n"
        "  --start-event   <num>    event number to start at when reading from a file (def: %d)\n"
        "  --input-file             file name for reading (def: %s)\n"
        "  --output-file            file name for writitng (def: %s)\n"
	"Combo spaghetti, that's with cole slaw: \n"
	"  --cmssw-simseeds         use CMS geom with simtracks for seeds \n"
	"  --cmssw-stdseeds         use CMS geom with CMSSW seeds uncleaned \n"
	"  --cmssw-n2seeds          use CMS geom with CMSSW seeds cleaned with N^2 routine \n"
	"  --cmssw-pureseeds        use CMS geom with pure CMSSW seeds (seeds which produced CMSSW reco tracks) \n"
	"  --cmssw-goodlabelseeds   use CMS geom with CMSSW seeds with label() >= 0 \n"
	"  --cmssw-val-trkparam     use CMSSW validation with track parameter matching \n"
	"  --cmssw-val-hit          use CMSSW validation with hit based matching (75 percent of reco track) \n"
	"  --cmssw-val-label        use CMSSW validation with track label matching (ONLY WITH PURE SEEDS!) \n"
        "GPU specific options: \n"
        "  --num-thr-reorg <num>    number of threads to run the hits reorganization (def: %d)\n"
        "  --seed-based             For CE. Switch to 1 CUDA thread per seed\n"
        "New options -- to be placed appropriately:\n"
        "  --kludge-cms-hit-errors  make sure err(xy) > 15 mum, err(z) > 30 mum (def: %s)\n"
        "  --backward-fit           perform backward fit during building (std and ce only) (def: %s)\n"
        ,
        argv[0],
        Config::geomPlugin.c_str(),
        Config::nEvents,
        Config::nTracks,
        Config::numThreadsSimulation, Config::numThreadsFinder, Config::numThreadsEvents,
        Config::chi2Cut,
	b2a(Config::usePhiQArrays),
        Config::numSeedsPerTask,
	Config::numHitsPerTask,
        Config::finderReportBestOutOfN,
	getOpt(Config::seedInput, g_seed_opts).c_str(),
	getOpt(Config::seedCleaning, g_clean_opts).c_str(),
        b2a(Config::cf_seeding),
        b2a(Config::cf_fitting),
	b2a(Config::readCmsswTracks),
        b2a(Config::readSimTrackStates),
        b2a(Config::quality_val),
        b2a(Config::root_val),
        b2a(Config::cmssw_val),
        b2a(Config::fit_val),
	b2a(Config::inclusiveShorts),
	getOpt(Config::cmsswMatching, g_match_opts).c_str(),
	b2a(Config::applyCMSSWHitMatch),
        b2a(Config::dumpForPlots),
        b2a(Config::silent),
        g_start_event,
      	g_input_file.c_str(),
      	g_output_file.c_str(),
	Config::numThreadsReorg,
        b2a(Config::kludgeCmsHitErrors),
        b2a(Config::backwardFit)
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
    else if (*i == "--use-phiq-arr")
    {
      Config::usePhiQArrays = true;
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
    else if(*i == "--build-fv")
    {
      g_run_build_all = false; g_run_build_fv = true;
    }
    else if (*i == "--seeds-per-task")
    {
      next_arg_or_die(mArgs, i);
      Config::numSeedsPerTask = atoi(i->c_str());
    }
    else if (*i == "--hits-per-task")
    {
      next_arg_or_die(mArgs, i);
      Config::numHitsPerTask = atoi(i->c_str());
    }
    else if(*i == "--best-out-of")
    {
      next_arg_or_die(mArgs, i);
      Config::finderReportBestOutOfN = atoi(i->c_str());
    }
    else if(*i == "--seed-input")
    {
      next_arg_or_die(mArgs, i);
      setOpt(*i,Config::seedInput,g_seed_opts,"seed input collection");
    }
    else if(*i == "--seed-cleaning")
    {
      next_arg_or_die(mArgs, i);
      setOpt(*i,Config::seedCleaning,g_clean_opts,"seed cleaning");
    }
    else if (*i == "--cf-seeding")
    {
      Config::cf_seeding = true;
    }
    else if (*i == "--cf-fitting")
    {
      Config::cf_fitting = true;
    }
    else if(*i == "--read-cmssw-tracks")
    {
      Config::readCmsswTracks = true;
    }
    else if (*i == "--read-simtrack-states")
    {
      Config::readSimTrackStates = true;
    }
    else if (*i == "--quality-val")
    {
      Config::quality_val = true; 
    }
    else if (*i == "--root-val")
    {
      Config::root_val = true; 
    }
    else if (*i == "--cmssw-val")
    {
      Config::cmssw_val = true;
    }
    else if (*i == "--fit-val")
    {
      Config::fit_val = true;
    }
    else if (*i == "--inc-shorts")
    {
      Config::inclusiveShorts = true;
    }
    else if(*i == "--cmssw-matching")
    {
      next_arg_or_die(mArgs, i);
      setOpt(*i,Config::cmsswMatching,g_match_opts,"CMSSW validation track matching");
    }
    else if (*i == "--hit-match")
    {
      Config::applyCMSSWHitMatch = true;
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
    else if (*i == "--start-event")
    {
      next_arg_or_die(mArgs, i);
      g_start_event = atoi(i->c_str());
    }
    else if (*i == "--input-file")
    {
      next_arg_or_die(mArgs, i);
      g_input_file = *i;
      g_operation = "read";
      Config::nEvents = -1;
    }
    else if (*i == "--output-file")
    {
      next_arg_or_die(mArgs, i);
      g_output_file = *i;
      g_operation = "write";
    }
    else if (*i == "--dump-for-plots")
    {
      Config::dumpForPlots = true;
    }
    else if (*i == "--silent")
    {
      Config::silent = true;
    }
    else if (*i == "--cmssw-simseeds")
    {
      Config::geomPlugin = "CMS-2017";
      Config::seedInput  = simSeeds;
    }
    else if (*i == "--cmssw-stdseeds")
    {
      Config::geomPlugin   = "CMS-2017";
      Config::seedInput    = cmsswSeeds;
      Config::seedCleaning = noCleaning;
    }
    else if (*i == "--cmssw-n2seeds")
    {
      Config::geomPlugin   = "CMS-2017";
      Config::seedInput    = cmsswSeeds;
      Config::seedCleaning = cleanSeedsN2;
    }
    else if (*i == "--cmssw-pureseeds")
    {
      Config::geomPlugin      = "CMS-2017";
      Config::seedInput       = cmsswSeeds;
      Config::seedCleaning    = cleanSeedsPure;
      Config::readCmsswTracks = true;
    }
    else if (*i == "--cmssw-goodlabelseeds")
    {
      Config::geomPlugin   = "CMS-2017";
      Config::seedInput    = cmsswSeeds;
      Config::seedCleaning = cleanSeedsBadLabel;
    }
    else if (*i == "--cmssw-val-trkparam")
    {
      Config::cmssw_val = true; 
      Config::readCmsswTracks = true;
      Config::cmsswMatching = trkParamBased;
    }
    else if (*i == "--cmssw-val-hit")
    {
      Config::cmssw_val = true; 
      Config::readCmsswTracks = true;
      Config::cmsswMatching = hitBased;
    }
    else if (*i == "--cmssw-val-label")
    {
      Config::cmssw_val = true; 
      Config::readCmsswTracks = true;
      Config::cmsswMatching = labelBased;
    }
    else if (*i == "--seed-based")
    {
      g_seed_based = true;
    }
    else if(*i == "--kludge-cms-hit-errors")
    {
      Config::kludgeCmsHitErrors = true;
    }
    else if(*i == "--backward-fit")
    {
      Config::backwardFit = true;
    }
    else
    {
      fprintf(stderr, "Error: Unknown option/argument '%s'.\n", i->c_str());
      exit(1);
    }

    mArgs.erase(start, ++i);
  }

  // Do some checking of options before going...
  if (Config::seedCleaning != cleanSeedsPure && Config::cmsswMatching == labelBased)
  {
    std::cerr << "What have you done?!? Can't mix cmssw label matching without pure seeds! Exiting..." << std::endl;
    exit(1);
  }

  // set to convert if I/O files both set!
  if (g_input_file != "" && g_output_file != "")
  {
    g_operation = "convert";
  }

  Config::RecalculateDependentConstants();

  printf ("Running with n_threads=%d, best_out_of=%d\n",
          Config::numThreadsFinder, Config::finderReportBestOutOfN);

  test_standard();

  return 0;
}
