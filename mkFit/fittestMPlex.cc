#include "fittest.h"

#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"
#include "Geometry.h"
//#define DEBUG
#include <Debug.h>

#include "MkFitter.h"
#if USE_CUDA
#include "fittestMPlex.h"
#include "FitterCU.h"
#include <omp.h>
#endif

#ifndef NO_ROOT
#include "TFile.h"
#include "TTree.h"
#include <mutex>
#endif

#include <tbb/tbb.h>

#include <iostream>
#include <memory>

#if defined(USE_VTUNE_PAUSE)
#include "ittnotify.h"
#endif

//==============================================================================
// runFittingTestPlex
//==============================================================================

#include "Pool.h"
namespace
{
  struct ExecutionContext
  {
    Pool<MkFitter>   m_fitters;

    void populate(int n_thr)
    {
      m_fitters.populate(n_thr - m_fitters.size());
    }
  };

  ExecutionContext g_exe_ctx;
  auto retfitr = [](MkFitter*   mkfp  ) { g_exe_ctx.m_fitters.ReturnToPool(mkfp);   };
}

double runFittingTestPlex(Event& ev, std::vector<Track>& rectracks)
{
   g_exe_ctx.populate(Config::numThreadsFinder);
   std::vector<Track>& simtracks = ev.simTracks_;

   const int Nhits = Config::nLayers;
   // XXX What if there's a missing / double layer?
   // Eventually, should sort track vector by number of hits!
   // And pass the number in on each "setup" call.
   // Reserves should be made for maximum possible number (but this is just
   // measurments errors, params).

   int theEnd = ( (Config::endcapTest && Config::readCmsswSeeds) ? ev.seedTracks_.size() : simtracks.size());
   int count = (theEnd + NN - 1)/NN;

#ifdef USE_VTUNE_PAUSE
   __itt_resume();
#endif

   double time = dtime();

   tbb::parallel_for(tbb::blocked_range<int>(0, count, std::max(1, Config::numSeedsPerTask/NN)),
     [&](const tbb::blocked_range<int>& i)
   {
     std::unique_ptr<MkFitter, decltype(retfitr)> mkfp(g_exe_ctx.m_fitters.GetFromPool(), retfitr);
     mkfp->SetNhits(Nhits);
     for (int it = i.begin(); it < i.end(); ++it)
     {
        int itrack = it * NN;
        int end    = itrack + NN;
	if (Config::endcapTest)
        {
	  //fixme, check usage of SlurpInTracksAndHits for endcapTest
	  if (Config::readCmsswSeeds) {
	    mkfp->InputSeedsTracksAndHits(ev.seedTracks_,simtracks, ev.layerHits_, itrack, end);
	  } else {
	    mkfp->InputTracksAndHits(simtracks, ev.layerHits_, itrack, end);
	  }
	  mkfp->FitTracksTestEndcap(end - itrack, &ev, true);
	}
        else
        {
        /*
         * MT, trying to slurp and fit at the same time ...
	  if (theEnd < end) {
	    end = theEnd;
	    mkfp->InputTracksAndHits(simtracks, ev.layerHits_, itrack, end);
	  } else {
	    mkfp->SlurpInTracksAndHits(simtracks, ev.layerHits_, itrack, end); // only safe for a full matriplex
	  }
	  
	  if (Config::cf_fitting) mkfp->ConformalFitTracks(true, itrack, end);
	  mkfp->FitTracks(end - itrack, &ev, true);
        */

          mkfp->InputTracksForFit(simtracks, itrack, end);

          // XXXX MT - for this need 3 points in ... right
	  // XXXX if (Config::cf_fitting) mkfp->ConformalFitTracks(true, itrack, end);

          mkfp->FitTracksWithInterSlurp(ev.layerHits_, end - itrack);
	}

	mkfp->OutputFittedTracks(rectracks, itrack, end);
     }
   });

   time = dtime() - time;

#ifdef USE_VTUNE_PAUSE
   __itt_pause();
#endif

   if (Config::fit_val) ev.Validate();

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
  separate_first_call_for_meaningful_profiling_numbers();

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
  int gplex_size = 10000;
  FitterCU<float> cuFitter(gplex_size);
  cuFitter.allocateDevice();
  cuFitter.allocate_extra_addBestHit();

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
        //if (g_run_fit_std) {
          std::string tree_name = "validation-plex-" + std::to_string(evt) + ".root";
        //}
      }
#endif
    }
    cuFitter.free_extra_addBestHit();
    cuFitter.freeDevice();
  }
  std::cerr << "###### [Fitting] Total GPU time: " << dtime() - total_gpu_time << " ######\n";
}


double runFittingTestPlexGPU(FitterCU<float> &cuFitter, 
    Event& ev, std::vector<Track>& rectracks)
{
  std::vector<Track>& simtracks = ev.simTracks_;

  cuFitter.createStream();

  Track *tracks_cu;
  cudaMalloc((void**)&tracks_cu, simtracks.size()*sizeof(Track));
  cudaMemcpyAsync(tracks_cu, &simtracks[0], simtracks.size()*sizeof(Track),
                  cudaMemcpyHostToDevice, cuFitter.get_stream());

  EventOfHitsCU events_of_hits_cu;
  events_of_hits_cu.allocGPU(ev.layerHits_);
  events_of_hits_cu.copyFromCPU(ev.layerHits_, cuFitter.get_stream());

  double time = dtime();

  cuFitter.FitTracks(tracks_cu, simtracks.size(), events_of_hits_cu, Config::nLayers);

  cudaMemcpy(&rectracks[0], tracks_cu, simtracks.size()*sizeof(Track), cudaMemcpyDeviceToHost);

  time = dtime() - time;


  events_of_hits_cu.deallocGPU();
  cudaFree(tracks_cu);

  cuFitter.destroyStream();

  return time;
}
#endif
