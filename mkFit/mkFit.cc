#include "MatriplexCommon.h"

#include "fittest.h"

#include "MkFitter.h"

#include "Timing.h"

#include <limits>

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

void test_standard()
{
  int  Ntracks  = 1024 * 1024;// * 1024; // * 10
  bool saveTree = false;

  generateTracks(simtracks, Ntracks);

  double tmp, tsm;

  smat_tracks.reserve(simtracks.size());
  tsm = runFittingTest(simtracks, smat_tracks);

  plex_tracks.resize(simtracks.size());
  tmp = runFittingTestPlex(simtracks, plex_tracks);

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

  printf("SMatrix = %.3f   Matriplex = %.3f   ---   SM/MP = %.3f\n", tsm, tmp, tsm / tmp);

#ifndef NO_ROOT
  make_validation_tree("validation-smat.root", simtracks, smat_tracks);
  make_validation_tree("validation-plex.root", simtracks, plex_tracks);
#endif
}

//==============================================================================

int main()
{
  // test_matriplex();

  // test_vtune();

  test_standard();

  return 0;
}
