#include "Matriplex/MatriplexCommon.h"

#include "fittestMPlex.h"
#include "buildtestMPlex.h"

#include "MkFitter.h"

#include "Config.h"

#include "Timing.h"

#include <limits>

#include "Event.h"

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
    float r = (l+1)*4.;
    VUSolid* utub = new VUSolid(r, r+.01);
    float z = r / std::tan(2.0*std::atan(std::exp(-eta))); // calculate z extent based on eta, r
    geom.AddLayer(utub, r, z);
  }
}

namespace
{
  FILE *g_file = 0;
  int   g_file_num_ev = 0;
  int   g_file_cur_ev = 0;

  std::string g_operation;
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

  for (int evt = 0; evt < Nevents; ++evt)
  {

    Event ev(geom, val, evt);

    ev.Simulate();//fixme add g_gen.seed(742331); and #pragma omp parallel for num_threads(NUM_THREADS_SIM)
    ev.resetLayerHitMap();

    fwrite(&Ntracks, sizeof(int), 1, fp);

    for (int i = 0; i < Ntracks; ++i)
    {
      ev.simTracks_[i].write_out(fp);
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
    simtracks[i].read_in(g_file);
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
  printf("  vusize=%i, num_th=%i\n",  MPT_SIZE, NUM_THREADS);
  printf("  sizeof(Track)=%zu, sizeof(Hit)=%zu, sizeof(SVector3)=%zu, sizeof(SMatrixSym33)=%zu, sizeof(MCHitInfo)=%zu\n",
         sizeof(Track), sizeof(Hit), sizeof(SVector3), sizeof(SMatrixSym33), sizeof(MCHitInfo));

  int Ntracks = Config::nTracks;
  // Ntracks  = 1;
  // bool saveTree = false;

  int Nevents = Config::nEvents;

  if (g_operation == "read")
  {
    Nevents = open_simtrack_file();
  }

  Geometry geom;
  initGeom(geom);
  Validation val;

  EventOfCombCandidates evocc;
  EventOfCombCandidates& event_of_comb_cands = evocc;

  EventOfCandidates evoc;
  EventOfCandidates& event_of_cands = evoc;

  double s_tmp=0, s_tsm=0, s_tsm2=0, s_tmp2=0, s_tsm2bh=0, s_tmp2bh=0;

  printf("\n");

  for (int evt = 1; evt <= Nevents; ++evt)
  {
    printf("Processing event %d\n", evt);

    Event ev(geom, val, evt);

    if (g_operation == "read")
    {
      Ntracks = read_simtrack_event(ev.simTracks_);
    }
    else
    {
      ev.Simulate();//fixme add g_gen.seed(742331); and #pragma omp parallel for num_threads(NUM_THREADS_SIM)
      ev.resetLayerHitMap();
    }

    plex_tracks.resize(ev.simTracks_.size());
    double tmp = runFittingTestPlex(ev, plex_tracks);

    double tmp2 = runBuildingTestPlex(ev, event_of_comb_cands);

    double tmp2bh = runBuildingTestPlexBestHit(ev, event_of_cands);

    printf("Matriplex fit = %.5f  --- Build  MX = %.5f  BHMX = %.5f\n",
           tmp, tmp2, tmp2bh);
    printf("\n");

    s_tmp    += tmp;
    s_tmp2   += tmp2;
    s_tmp2bh += tmp2bh;
  }
  printf("================================================================\n");
  printf("=== TOTAL for %d events\n", Nevents);
  printf("================================================================\n");

  printf("Matriplex fit = %.5f  --- Build  MX = %.5f  BHMX = %.5f\n",
         s_tmp, s_tmp2, s_tmp2bh);

  if (g_operation == "read")
  {
    close_simtrack_file();
  }

#ifndef NO_ROOT
  make_validation_tree("validation-plex.root", ev.simTracks_, plex_tracks);
#endif
}

//==============================================================================

void usage_and_die(const char* name)
{
  fprintf(stderr,
          "Usage:\n"
          "  %s                  --> runs simulation between events\n"
          "  %s write [filename] --> runs simulation only, outputs events to file\n"
          "  %s read  [filename] --> runs reco only, reads events from file\n"
          "Default filename is \"simtracks.bin\".\n", name, name, name);
  exit(1);
}


int main(int argc, const char *argv[])
{
#ifdef USE_VTUNE_PAUSE
  __itt_pause();
#endif

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

  // test_matriplex();

  // test_vtune();

  if (g_operation == "write")
  {
    generate_and_save_tracks();
  }
  else
  {
    test_standard();
  }

  return 0;
}
