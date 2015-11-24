#include "fittest.h"

#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"
#include "Geometry.h"

#include "MkFitter.h"

#ifndef NO_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

#include <omp.h>

#include <iostream>
#include <memory>

//==============================================================================

void make_validation_tree(const char         *fname,
                          std::vector<Track> &simtracks,
                          std::vector<Track> &rectracks)
{
#ifndef NO_ROOT

   assert(simtracks.size() == rectracks.size());

   float pt_mc, pt_fit, pt_err, chg;

   TFile *file = TFile::Open(fname, "recreate");
   TTree *tree = new TTree("T", "Validation Tree for simple Kalman fitter");;

   tree->Branch("pt_mc",  &pt_mc,  "pt_mc");
   tree->Branch("pt_fit", &pt_fit, "pt_fit");
   tree->Branch("pt_err", &pt_err, "pt_err");
   tree->Branch("chg",    &chg,    "chg");

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

      tree->Fill();
   }

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

   double time = dtime();

#pragma omp parallel for
   for (int itrack = 0; itrack < theEnd; itrack += NN)
   {
      int end = std::min(itrack + NN, theEnd);

      MkFitter *mkfp = mkfp_arr[omp_get_thread_num()];

      mkfp->InputTracksAndHits(simtracks, ev.layerHits_, itrack, end);

      mkfp->FitTracks();

#ifndef NO_ROOT
      mkfp->OutputFittedTracks(rectracks, itrack, end);
#endif
   }

   time = dtime() - time;

   for (int i = 0; i < Config::numThreadsFinder; ++i)
   {
     _mm_free(mkfp_arr[i]);
   }
   //_mm_free(mkfp);

   return time;
}
