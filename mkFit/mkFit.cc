#include "Matriplex/MatriplexCommon.h"

#include "fittestMPlex.h"
#include "buildtestMPlex.h"

#include "MkFitter.h"

#include "Timing.h"

#include <limits>

#if defined(USE_VTUNE_PAUSE)
#include "ittnotify.h"
#endif

std::vector<Track> simtracks;

std::vector<Track> smat_tracks;
std::vector<Track> plex_tracks;

//==============================================================================

MkFitter *g_mkfp;

const int Nhits = MAX_HITS; // XXXXX ARGH !!!! What if there's a missing / double layer?

const int Nloop = 100;

//==============================================================================

long64 single_run(int                 n_tracks,
                  MkFitter           *mkfp,
                  std::vector<Track> &trk_in,
                  std::vector<Track> &trk_out)
{
  int theEnd = n_tracks;

  for (int itrack = 0; itrack < theEnd; itrack += NN)
  {
    int end = std::min(itrack + NN, theEnd);

    mkfp->InputTracksAndHits(trk_in, itrack, end);

    // Store input track params into a stash for faster reinitialization.
    // Charge doesn't get changed.
    // If you activate this, also swap the reinit in the loop below.
    // This doesn't really help all that much, it seems.
    //
    // MPlexLS err; err = mkfp->GetErr0();
    // MPlexLV par; par = mkfp->GetPar0();

    for (int x = 0; x < Nloop; ++x)
    {
      if (x != 0)
      {
        mkfp->InputTracksOnly(trk_in, itrack, end);
        // mkfp->GetErr0() = err;
        // mkfp->GetPar0() = par;
      }

      mkfp->FitTracks();
    }

#ifndef NO_ROOT
    mkfp->OutputFittedTracks(trk_out, itrack, end);
#endif
  }

  // For propagateLine
  // return long64(Nloop) * n_tracks * (68 + 306) * Nhits;

  return long64(Nloop) * n_tracks * (1200 + 306) * Nhits;
}

//------------------------------------------------------------------------------

long64 single_run_glob(int n_tracks)
{
  return single_run(n_tracks, g_mkfp, simtracks, plex_tracks);
}

//==============================================================================

void test_matriplex()
{
  int Nmin  = 16;
  int Nmax  = 64 * 1024; // 32 * 1024;

  generateTracks(simtracks, Nmax);

  g_mkfp = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);

  g_mkfp->CheckAlignment();

  Timing t([&](int n_vec)
           {
             return single_run_glob(n_vec);
           });

  t.print_tuple_header();

  // Warm up
  t.time_loop(512, 4);

  for (int n = Nmin; n <= Nmax; n *= 2)
  {
    // Nloop = std::numeric_limits<int>::max() / (374 * Nhits * n);

    // printf("XXX n=%d, nloop=%d\n", n, Nloop);

    t.auto_time_loop(n, 2);
    
    // t.time_loop(n, 1);

    t.print(n);
  }

  _mm_free(g_mkfp);
}

//==============================================================================

void test_vtune()
{
  int Nmax  = 512;
  // 512 here is calculated so the whole thing fits into L2 cache.
  // For single thread per core one could go to twice that.

  // KMP_AFFINITY=scatter KMP_PLACE_THREADS=1T
  // KMP_AFFINITY=compact KMP_PLACE_THREADS=2T

  // #pragma omp parallel for num_threads(NUM_THREADS)
  // for (int i = 0; i < NUM_THREADS; ++i)
  {
    std::vector<Track> sim_trk;
    std::vector<Track> rec_trk;

    generateTracks(sim_trk, Nmax);
    rec_trk.resize(Nmax);

    MkFitter *mf = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);

    mf->CheckAlignment();

    for (int i = 0 ; i < 200; ++i) 
    {
      single_run(Nmax, mf, sim_trk, rec_trk);
    }

    _mm_free(mf);
  }
}

//==============================================================================

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

  int Ntracks = Config::g_NTracks;

  int Nevents = Config::g_NEvents;

  fwrite(&Nevents, sizeof(int), 1, fp);

  for (int ev = 0; ev < Nevents; ++ev)
  {
    generateTracks(simtracks, Ntracks);

    fwrite(&Ntracks, sizeof(int), 1, fp);

    for (int i = 0; i < Ntracks; ++i)
    {
      simtracks[i].write_out(fp);
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
  //   nb = Config::getBothEtaBins(eta, b1, b2);
  //   printf("eta=%6.2f  bin=%3d  bin1=%3d  bin2=%3d nb=%d\n",
  //          eta, Config::getEtaBin(eta), b1, b2, nb);
  // }

  // return;
  // ---- end MT test

  int Ntracks = Config::g_NTracks;
  // Ntracks  = 1;
  // bool saveTree = false;

  int Nevents = Config::g_NEvents;

  if (g_operation == "read")
  {
    Nevents = open_simtrack_file();
  }

  double s_tmp=0, s_tsm=0, s_tsm2=0, s_tmp2=0, s_tsm2bh=0, s_tmp2bh=0;

  for (int ev = 0; ev < Nevents; ++ev)
  {
    if (g_operation == "read")
    {
      Ntracks = read_simtrack_event(simtracks);
    }
    else
    {
      generateTracks(simtracks, Ntracks);
    }

    double tmp, tsm;

    smat_tracks.reserve(simtracks.size());
    tsm = 0; // runFittingTest(simtracks, smat_tracks);

    plex_tracks.resize(simtracks.size());
    tmp = 0; // runFittingTestPlex(simtracks, plex_tracks);

    double tsm2 = 0;//runBuildingTest(simtracks);
    double tmp2 = runBuildingTestPlex(simtracks);

    double tsm2bh = 0;//runBuildingTestBestHit(simtracks);
    double tmp2bh = 0;//runBuildingTestPlexBestHit(simtracks);

    // Second pass -- select problematic tracks and refit them
    if (false)
    {
      int iout = 0;
      for (int i = 0; i < Ntracks; ++i)
      {
        SVector6 &simp = simtracks[i].parameters();
        SVector6 &recp = plex_tracks[i].parameters();

        float pt_mc  = sqrt(simp[3]*simp[3] + simp[4]*simp[4]);
        float pt_fit = sqrt(recp[3]*recp[3] + recp[4]*recp[4]);

        if (std::abs((pt_mc - pt_fit) / pt_mc) > 100)
        {
          printf("Got bad track: %d %d %f %f\n", i, iout, pt_mc, pt_fit);
          if (i != 0)
            simtracks[iout] = simtracks[i];
          ++iout;
          if (iout >= 16)
            break;
        }
      }

      g_dump = true;

      simtracks.resize(16);
      smat_tracks.resize(0); smat_tracks.reserve(16);
      plex_tracks.resize(16);

      tsm = runFittingTest(simtracks, smat_tracks);

      printf("\n\n\n===========================================================\n\n\n");

      tmp = runFittingTestPlex(simtracks, plex_tracks);
    }

    printf("SMatrix = %.5f   Matriplex = %.5f   ---   SM/MP = %.5f  --- Build SM = %.5f    MX = %.5f    BHSM = %.5f    BHMX = %.5f\n",
           tsm, tmp, tsm / tmp, tsm2, tmp2, tsm2bh, tmp2bh);

    s_tmp    += tmp;    s_tsm    += tsm;
    s_tsm2   += tsm2;   s_tmp2   += tmp2;
    s_tsm2bh += tsm2bh; s_tmp2bh += tmp2bh;
  }
  printf("================================================================\n");
  printf("=== TOTAL for %d events\n", Nevents);
  printf("================================================================\n");

  printf("SMatrix = %.5f   Matriplex = %.5f   ---   SM/MP = %.5f  --- Build SM = %.5f    MX = %.5f    BHSM = %.5f    BHMX = %.5f\n",
         s_tsm, s_tmp, s_tsm / s_tmp, s_tsm2, s_tmp2, s_tsm2bh, s_tmp2bh);

  if (g_operation == "read")
  {
    close_simtrack_file();
  }

#ifndef NO_ROOT
  make_validation_tree("validation-smat.root", simtracks, smat_tracks);
  make_validation_tree("validation-plex.root", simtracks, plex_tracks);
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
