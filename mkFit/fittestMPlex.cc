#include "fittest.h"

#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"
#include "Geometry.h"

#include "MkFitter.h"
#if USE_CUDA
#include "fittestMPlex.h"
#include "FitterCU.h"
#endif

#ifndef NO_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

//#define DEBUG
#include <Debug.h>

#include <omp.h>

#include <iostream>
#include <memory>
#include <mutex>

#if defined(USE_VTUNE_PAUSE)
#include "ittnotify.h"
#endif

//==============================================================================

void make_validation_tree(const char         *fname,
                          std::vector<Track> &simtracks,
                          std::vector<Track> &rectracks)
{
   assert(simtracks.size() == rectracks.size());

   float pt_mc, pt_fit, pt_err, chg;
   int goodtrk = 0;

#ifndef NO_ROOT
   static std::mutex roolock;

   std::lock_guard<std::mutex> rooguard(roolock);

   TFile *file = TFile::Open(fname, "recreate");
   TTree *tree = new TTree("T", "Validation Tree for simple Kalman fitter");;

   tree->Branch("pt_mc",  &pt_mc,  "pt_mc");
   tree->Branch("pt_fit", &pt_fit, "pt_fit");
   tree->Branch("pt_err", &pt_err, "pt_err");
   tree->Branch("chg",    &chg,    "chg");
#endif

#ifdef USE_CUDA
   std::vector<float> diff_vec;
#endif

   const int NT = simtracks.size();
   for (int i = 0; i < NT; ++i)
   {
      SVector6     &simp   = simtracks[i].parameters_nc();
      SVector6     &recp   = rectracks[i].parameters_nc();
      SMatrixSym66 &recerr = rectracks[i].errors_nc();

      pt_mc  = sqrt(simp[3]*simp[3] + simp[4]*simp[4]);
      pt_fit = sqrt(recp[3]*recp[3] + recp[4]*recp[4]);
      pt_err = sqrt(recerr[3][3]*recp[3]*recp[3] +
                    recerr[4][4]*recp[4]*recp[4] +
                    recerr[3][4]*recp[3]*recp[4] * 2) / pt_fit;
      chg = simtracks[i].charge();


#ifndef NO_ROOT
      tree->Fill();
#endif

#ifdef USE_CUDA
      float diff = (pt_mc - pt_fit) / pt_err;
      if (std::abs(diff) < 5.0f) {
        diff_vec.push_back(diff);
        ++goodtrk;
      } else {
        dprint("pt_mc, pt_fit, pt_err " << pt_mc << " " << pt_fit << " " << pt_err);
      }
#endif
   }
#ifdef USE_CUDA
   float mean = std::accumulate(diff_vec.begin(), diff_vec.end(), 0.0)
              / diff_vec.size();

   std::transform(diff_vec.begin(), diff_vec.end(), 
                  diff_vec.begin(),  // in-place
                  [mean](float x) {return (x-mean)*(x-mean);});
                  
   float stdev = std::sqrt(
       std::accumulate(diff_vec.begin(), diff_vec.end(), 0.0)
       / diff_vec.size());

   std::cerr << goodtrk << " good tracks, mean pt pull: " << mean
             << " standard dev: " << stdev << std::endl;
#endif

#ifndef NO_ROOT
   file->Write();
   file->Close();
   delete file;
#endif
}

//==============================================================================
// runFittingTestPlex
//==============================================================================

double runFittingTestPlex(Event& ev, std::vector<Track>& rectracks)
{

   std::vector<Track>& simtracks = ev.simTracks_;

   const int Nhits = MAX_HITS;
   // XXX What if there's a missing / double layer?
   // Eventually, should sort track vector by number of hits!
   // And pass the number in on each "setup" call.
   // Reserves should be made for maximum possible number (but this is just
   // measurments errors, params).

   // NOTE: MkFitter *MUST* be on heap, not on stack!
   // Standard operator new screws up alignment of ALL MPlex memebrs of MkFitter,
   // even if one adds attr(aligned(64)) thingy to every possible place.

   // MkFitter *mkfp = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);

   std::vector<MkFitter*> mkfp_arr(Config::numThreadsFinder);

   for (int i = 0; i < Config::numThreadsFinder; ++i)
   {
     mkfp_arr[i] = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);
   }

   int theEnd = simtracks.size();

#ifdef USE_VTUNE_PAUSE
   __itt_resume();
#endif

   double time = dtime();

#pragma omp parallel for
   for (int itrack = 0; itrack < theEnd; itrack += NN)
   {
      int end = std::min(itrack + NN, theEnd);

      MkFitter *mkfp = mkfp_arr[omp_get_thread_num()];

      //mkfp->InputTracksAndHits(simtracks, ev.layerHits_, itrack, end);
      mkfp->SlurpInTracksAndHits(simtracks, ev.layerHits_, itrack, end);

      if (Config::cf_fitting) mkfp->ConformalFitTracks(true, itrack, end);
      mkfp->FitTracks(end - itrack);

#ifndef NO_ROOT
      mkfp->OutputFittedTracks(rectracks, itrack, end);
#endif
   }

   time = dtime() - time;

#ifdef USE_VTUNE_PAUSE
   __itt_pause();
#endif

   for (int i = 0; i < Config::numThreadsFinder; ++i)
   {
     _mm_free(mkfp_arr[i]);
   }
   //_mm_free(mkfp);

   return time;
}

#ifdef USE_CUDA
void runAllEventsFittingTestPlexGPU(std::vector<Event>& events)
{
  double s_tmp = 0.0;
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

      //if (g_run_fit_std) tmp = runFittingTestPlexGPU(cuFitter, ev, plex_tracks_ev);
      runFittingTestPlexGPU(cuFitter, ev, plex_tracks_ev);

      printf("Matriplex fit = %.5f  -------------------------------------", tmp);
      printf("\n");
      s_tmp    += tmp;
#if 0  // 0 for timing, 1 for validation
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
  std::cerr << "###### [Fitting] Total GPU time: " << dtime() - total_gpu_time << " ######\n";
}


double runFittingTestPlexGPU(FitterCU<float> &cuFitter, 
    Event& ev, std::vector<Track>& rectracks)
{

   std::vector<Track>& simtracks = ev.simTracks_;

   const int Nhits = MAX_HITS;
   // XXX What if there's a missing / double layer?
   // Eventually, should sort track vector by number of hits!
   // And pass the number in on each "setup" call.
   // Reserves should be made for maximum possible number (but this is just
   // measurments errors, params).

   // NOTE: MkFitter *MUST* be on heap, not on stack!
   // Standard operator new screws up alignment of ALL MPlex memebrs of MkFitter,
   // even if one adds attr(aligned(64)) thingy to every possible place.

   // MkFitter *mkfp = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);

   MkFitter* mkfp_arr = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);

   int theEnd = simtracks.size();
   double time = dtime();
   int Nstride = NN;

   for (int itrack = 0; itrack < theEnd; itrack += Nstride)
   {
      int end = std::min(itrack + Nstride, theEnd);

      MkFitter *mkfp = mkfp_arr;

      //double time_input = dtime();
      mkfp->InputTracksAndHits(simtracks, ev.layerHits_, itrack, end);
      //std::cerr << "Input time: " << (dtime() - time_input)*1e3 << std::endl;

      cuFitter.FitTracks(mkfp->Chg,
                         mkfp->GetPar0(),
                         mkfp->GetErr0(),
                         mkfp->msPar,
                         mkfp->msErr,
                         Nhits,
                         simtracks, itrack, end, ev.layerHits_);

//#ifndef NO_ROOT
      double time_output = dtime();
      mkfp->OutputFittedTracks(rectracks, itrack, end);
      std::cerr << "Output time: " << (dtime() - time_output)*1e3 << std::endl;
//#endif
   }

   time = dtime() - time;

   _mm_free(mkfp_arr);

   return time;
}
#endif
