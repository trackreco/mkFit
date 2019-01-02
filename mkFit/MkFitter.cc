#include "MkFitter.h"

#include "KalmanUtilsMPlex.h"
#include "ConformalUtilsMPlex.h"
#include "MatriplexPackers.h"

#ifdef USE_CUDA
//#include "FitterCU.h"
#endif

//#define DEBUG
#include "Debug.h"

#include <sstream>

namespace mkfit {

void MkFitter::CheckAlignment()
{
  printf("MkFitter alignment check:\n");
  Matriplex::align_check("  Err[0]   =", &Err[0].fArray[0]);
  Matriplex::align_check("  Err[1]   =", &Err[1].fArray[0]);
  Matriplex::align_check("  Par[0]   =", &Par[0].fArray[0]);
  Matriplex::align_check("  Par[1]   =", &Par[1].fArray[0]);
  Matriplex::align_check("  msErr[0] =", &msErr[0].fArray[0]);
  Matriplex::align_check("  msPar[0] =", &msPar[0].fArray[0]);
}

void MkFitter::PrintPt(int idx)
{
  for (int i = 0; i < NN; ++i)
  {
    printf("%5.2f  ", std::hypot(Par[idx].At(i, 3, 0), Par[idx].At(i, 4, 0)));
  }
}

int MkFitter::countValidHits(int itrack, int end_hit) const
{
  int result = 0;
  for (int hi = 0; hi < end_hit; ++hi)
    {
      if (HoTArr[hi](itrack, 0, 0).index >= 0) result++;
    }
  return result;
}

int MkFitter::countInvalidHits(int itrack, int end_hit) const
{
  int result = 0;
  for (int hi = 0; hi < end_hit; ++hi)
    {
      // XXXX MT: Should also count -2 hits as invalid?
      if (HoTArr[hi](itrack, 0, 0).index == -1) result++;
    }
  return result;
}

//==============================================================================

void MkFitter::InputTracksAndHits(const std::vector<Track>&  tracks,
                                  const std::vector<HitVec>& layerHits,
                                  int beg, int end)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;

// FIXME: uncomment when track building is ported to GPU.
#if USE_CUDA_NOT_YET
//#ifdef USE_CUDA
  // This openmp loop brings some performances when using
  // a single thread to fit all events.
  // However, it is more advantageous to use the threads to
  // parallelize over Events.
  omp_set_num_threads(Config::numThreadsReorg);
#pragma omp parallel for private(itrack)
#endif
  for (int i = beg; i < end; ++i, ++itrack)
  {
    const Track &trk = tracks[i];

    Err[iC].CopyIn(itrack, trk.errors().Array());
    Par[iC].CopyIn(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0)   = trk.charge();
    Chi2(itrack, 0, 0)  = trk.chi2();
    Label(itrack, 0, 0) = trk.label();

// CopyIn seems fast enough, but indirections are quite slow.
// For GPU computations, it has been moved in between kernels
// in an attempt to overlap CPU and GPU computations.
// FIXME: uncomment when track building is ported to GPU.
#if 1
//#ifndef USE_CUDA
    for (int hi = 0; hi < Nhits; ++hi)
    {
      HoTArr[hi](itrack, 0, 0) = trk.getHitOnTrack(hi);

      const int hidx = trk.getHitIdx(hi);
      if (hidx < 0) continue;

      const Hit &hit = layerHits[hi][hidx];
      msErr[hi].CopyIn(itrack, hit.errArray());
      msPar[hi].CopyIn(itrack, hit.posArray());
    }
#endif
  }
}

void MkFitter::InputTracksAndHits(const std::vector<Track>&  tracks,
                                  const std::vector<LayerOfHits>& layerHits,
                                  int beg, int end)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack;
//#ifdef USE_CUDA
#if 0
  // This openmp loop brings some performances when using
  // a single thread to fit all events.
  // However, it is more advantageous to use the threads to
  // parallelize over Events.
  omp_set_num_threads(Config::numThreadsReorg);
#pragma omp parallel for private(itrack)
#endif
  for (int i = beg; i < end; ++i) {
    itrack = i - beg;
    const Track &trk = tracks[i];

    Label(itrack, 0, 0) = trk.label();

    Err[iC].CopyIn(itrack, trk.errors().Array());
    Par[iC].CopyIn(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();

// CopyIn seems fast enough, but indirections are quite slow.
// For GPU computations, it has been moved in between kernels
// in an attempt to overlap CPU and GPU computations.
//#ifndef USE_CUDA
#if 1
    for (int hi = 0; hi < Nhits; ++hi)
    {
      const int hidx = trk.getHitIdx(hi);
      const int hlyr = trk.getHitLyr(hi);
      const Hit &hit = layerHits[hlyr].m_hits[hidx];

      msErr[hi].CopyIn(itrack, hit.errArray());
      msPar[hi].CopyIn(itrack, hit.posArray());

      HoTArr[hi](itrack, 0, 0) = trk.getHitOnTrack(hi);
    }
#endif
  }
}

void MkFitter::SlurpInTracksAndHits(const std::vector<Track>&  tracks,
                                    const std::vector<HitVec>& layerHits,
                                    int beg, int end)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  MatriplexTrackPacker mtp(tracks[beg]);

//#ifdef USE_CUDA
#if 0
  // This openmp loop brings some performances when using
  // a single thread to fit all events.
  // However, it is more advantageous to use the threads to
  // parallelize over Events.
  omp_set_num_threads(Config::numThreadsReorg);
#pragma omp parallel for private(itrack)
#endif
  for (int i = beg; i < end; ++i)
  {
    int itrack = i - beg;
    const Track &trk = tracks[i];

    Label(itrack, 0, 0) = trk.label();

    mtp.AddInput(trk);

    Chg(itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();
  }

  mtp.Pack(Err[iC], Par[iC]);

// CopyIn seems fast enough, but indirections are quite slow.
// For GPU computations, it has been moved in between kernels
// in an attempt to overlap CPU and GPU computations.
//#ifndef USE_CUDA
#if 1

  for (int hi = 0; hi < Nhits; ++hi)
  {
    MatriplexHitPacker mhp(layerHits[hi][0]);

    for (int i = beg; i < end; ++i)
    {
      const int   hidx = tracks[i].getHitIdx(hi);
      const Hit  &hit  = layerHits[hi][hidx];

      HoTArr[hi](i - beg, 0, 0) = tracks[i].getHitOnTrack(hi);

      mhp.AddInput(hit);
    }

    mhp.Pack(msErr[hi], msPar[hi]);
  }
#endif
}

void MkFitter::InputTracksAndHitIdx(const std::vector<Track>& tracks,
                                    int beg, int end,
                                    bool inputProp)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  const int iI = inputProp ? iP : iC;

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    const Track &trk = tracks[i];

    Err[iI].CopyIn(itrack, trk.errors().Array());
    Par[iI].CopyIn(itrack, trk.parameters().Array());

    Chg (itrack, 0, 0)  = trk.charge();
    Chi2(itrack, 0, 0)  = trk.chi2();
    Label(itrack, 0, 0) = trk.label();

    for (int hi = 0; hi < Nhits; ++hi)
    {
      HoTArr[hi](itrack, 0, 0) = trk.getHitOnTrack(hi);
    }
  }
}

void MkFitter::InputTracksAndHitIdx(const std::vector<std::vector<Track> >& tracks,
                                    const std::vector<std::pair<int,int> >& idxs,
                                    int beg, int end, bool inputProp)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  const int iI = inputProp ? iP : iC;

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    const Track &trk = tracks[idxs[i].first][idxs[i].second];

    Label(itrack, 0, 0) = trk.label();
    SeedIdx(itrack, 0, 0) = idxs[i].first;
    CandIdx(itrack, 0, 0) = idxs[i].second;

    Err[iI].CopyIn(itrack, trk.errors().Array());
    Par[iI].CopyIn(itrack, trk.parameters().Array());

    Chg (itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();

    for (int hi = 0; hi < Nhits; ++hi)
    {
      HoTArr[hi](itrack, 0, 0) = trk.getHitOnTrack(hi);
    }
  }
}

void MkFitter::InputSeedsTracksAndHits(const std::vector<Track>&  seeds,
				       const std::vector<Track>&  tracks,
				       const std::vector<HitVec>& layerHits,
				       int beg, int end)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack;
  for (int i = beg; i < end; ++i)
  {
    itrack = i - beg;

    const Track &see = seeds[i];

    Label(itrack, 0, 0) = see.label();
    if (see.label()<0) continue;

    Err[iC].CopyIn(itrack, see.errors().Array());
    Par[iC].CopyIn(itrack, see.parameters().Array());

    Chg(itrack, 0, 0) = see.charge();
    Chi2(itrack, 0, 0) = see.chi2();

    const Track &trk = tracks[see.label()];

// CopyIn seems fast enough, but indirections are quite slow.
// For GPU computations, it has been moved in between kernels
// in an attempt to overlap CPU and GPU computations.
//#ifndef USE_CUDA
#if 1
    for (int hi = 0; hi < Nhits; ++hi)
    {
      HoTArr[hi](itrack, 0, 0) = trk.getHitOnTrack(hi);

      const int hidx = trk.getHitIdx(hi);
      if (hidx < 0) continue; //fixme, check if this is harmless

      const Hit &hit = layerHits[hi][hidx];
      msErr[hi].CopyIn(itrack, hit.errArray());
      msPar[hi].CopyIn(itrack, hit.posArray());
    }
#endif
  }
}

//------------------------------------------------------------------------------
// Fitting with interleaved hit loading
//------------------------------------------------------------------------------

void MkFitter::InputTracksForFit(const std::vector<Track>& tracks,
                                 int beg, int end)
{
  // Loads track parameters and hit indices.

  // XXXXMT4K has Config::nLayers: How many hits do we read in?
  // Check for max? Expect an argument?
  // What to do with invalid hits? Skip?

  // XXXX MT Here the same idx array WAS used for slurping in of tracks and
  // hots. With this, two index arrays are built, one within each packer.

  MatriplexTrackPacker mtp  ( tracks[beg]);
  MatriplexHoTPacker   mhotp(*tracks[beg].getHitsOnTrackArray());

  int itrack = 0;

  for (int i = beg; i < end; ++i, ++itrack)
  {
    const Track &trk = tracks[i];

    Chg(itrack, 0, 0)   = trk.charge();
    Chi2(itrack, 0, 0)  = trk.chi2();
    Label(itrack, 0, 0) = trk.label();

    mtp.AddInput(trk);

    mhotp.AddInput( * trk.getHitsOnTrackArray());
  }

  mtp.Pack(Err[iC], Par[iC]);
  for (int ll = 0; ll < Config::nLayers; ++ll)
  {
     mhotp.Pack(HoTArr[ll], ll);
  }
}

void MkFitter::FitTracksWithInterSlurp(const std::vector<HitVec>& layersohits,
                                       const int N_proc)
{
  // XXXX This has potential issues hits coming from different layers!
  // Expected to only work reliably with barrel (consecutive layers from 0 -> Nhits)
  // and with hits present on every layer for every track.

  // Loops over layers and:
  // a) slurps in hit parameters;
  // b) propagates and updates tracks

  for (int ii = 0; ii < Nhits; ++ii)
  {
    // XXXX Assuming hit index corresponds to layer!
    MatriplexHitPacker mhp(layersohits[ii][0]);

    for (int i = 0; i < N_proc; ++i)
    {
      const int hidx = HoTArr[ii](i, 0, 0).index;
      const int hlyr = HoTArr[ii](i, 0, 0).layer;

      // XXXXMT4K What to do with hidx < 0 ????
      // This could solve the unbalanced fit.
      // Or, if the hidx is the "universal" missing hit, it could just work.
      // Say, hidx = 0 ... grr ... but then we don't know it is missing.

      if (hidx < 0 || hlyr < 0) {
         mhp.AddNullInput();
      } else {
         mhp.AddInput( layersohits[hlyr][hidx] );
      }
    }

    mhp.Pack(msErr[0], msPar[0]);

    PropagateTracksToHitR(msPar[0], N_proc, Config::forward_fit_pflags);

    kalmanUpdate(Err[iP], Par[iP], msErr[0], msPar[0],
                 Err[iC], Par[iC], N_proc);
  }
}

//==============================================================================
// Fitting functions
//==============================================================================

void MkFitter::ConformalFitTracks(bool fitting, int beg, int end)
{
  // bool fitting to determine to use fitting CF error widths
  // in reality, this is depedent on hits used to make pulls 
  // could consider writing an array for widths for a given hit combo 
  // to give precise widths --> then would drop boolean
  // also used to determine which hits to use

  int front,middle,back;
  
  // FIXME FITTING HITS --> assume one hit per layer and all layers found! BAD! Need vector of indices to do this right instead... 
  // can always assume 0,1,2 for seeding --> triplets in forward direction
#ifdef INWARDFIT
  front  = (fitting ?  Config::nLayers-1    : 0); // i.e. would rather have true option not hardcoded... but set by ACTUAL last hit found
  middle = (fitting ? (Config::nLayers-1)/2 : 1); // same with this one... would rather middle hit be in the middle!
  back   = (fitting ? 0 : 2); 
#else
  front  = (fitting ? 0 : 0); 
  middle = (fitting ? (Config::nLayers-1)/2 : 1); // ditto above
  back   = (fitting ?  Config::nLayers-1    : 2); // yup...
#endif

  // write to iC --> next step will be a propagation no matter what
  conformalFitMPlex(fitting, Label, Err[iC], Par[iC], 
		    msPar[front], msPar[middle], msPar[back]);

  // need to set most off-diagonal elements in unc. to zero, inflate all other elements;
  if (fitting) 
  { 
  using idx_t = Matriplex::idx_t;
  const idx_t N = NN;
#pragma omp simd
    for (int n = 0; n < N; ++n)
    {
      Err[iC].At(n, 0, 0) = Err[iC].ConstAt(n, 0, 0) * Config::blowupfit; 
      Err[iC].At(n, 0, 1) = Err[iC].ConstAt(n, 0, 1) * Config::blowupfit; 
      Err[iC].At(n, 1, 0) = Err[iC].ConstAt(n, 1, 0) * Config::blowupfit; 
      Err[iC].At(n, 1, 1) = Err[iC].ConstAt(n, 1, 1) * Config::blowupfit; 
      Err[iC].At(n, 2, 2) = Err[iC].ConstAt(n, 2, 2) * Config::blowupfit;
      Err[iC].At(n, 3, 3) = Err[iC].ConstAt(n, 3, 3) * Config::blowupfit;
      Err[iC].At(n, 4, 4) = Err[iC].ConstAt(n, 4, 4) * Config::blowupfit;
      Err[iC].At(n, 5, 5) = Err[iC].ConstAt(n, 5, 5) * Config::blowupfit;
      
      Err[iC].At(n, 0, 2) = 0.0f; Err[iC].At(n, 0, 3) = 0.0f; 
      Err[iC].At(n, 0, 4) = 0.0f; Err[iC].At(n, 0, 5) = 0.0f;
      Err[iC].At(n, 1, 2) = 0.0f; Err[iC].At(n, 1, 3) = 0.0f; 
      Err[iC].At(n, 1, 4) = 0.0f; Err[iC].At(n, 1, 5) = 0.0f;      
      Err[iC].At(n, 2, 0) = 0.0f; Err[iC].At(n, 2, 1) = 0.0f;
      Err[iC].At(n, 2, 3) = 0.0f; Err[iC].At(n, 2, 4) = 0.0f; 
      Err[iC].At(n, 2, 5) = 0.0f; Err[iC].At(n, 3, 0) = 0.0f; 
      Err[iC].At(n, 3, 1) = 0.0f; Err[iC].At(n, 3, 2) = 0.0f;
      Err[iC].At(n, 3, 4) = 0.0f; Err[iC].At(n, 3, 5) = 0.0f;
      Err[iC].At(n, 4, 0) = 0.0f; Err[iC].At(n, 4, 1) = 0.0f; 
      Err[iC].At(n, 4, 2) = 0.0f; Err[iC].At(n, 4, 3) = 0.0f; 
      Err[iC].At(n, 4, 5) = 0.0f; Err[iC].At(n, 5, 0) = 0.0f; 
      Err[iC].At(n, 5, 1) = 0.0f; Err[iC].At(n, 5, 2) = 0.0f;
      Err[iC].At(n, 5, 3) = 0.0f; Err[iC].At(n, 5, 4) = 0.0f;
    }
  }
}

void MkFitter::FitTracks(const int N_proc, const Event * ev, const PropagationFlags pflags)
{
  // Fitting loop.

  for (int hi = 0; hi < Nhits; ++hi)
  {
    // Note, charge is not passed (line propagation).
    // propagateLineToRMPlex(Err[iC], Par[iC], msErr[hi], msPar[hi],
    //                       Err[iP], Par[iP]);

    PropagateTracksToHitR(msPar[hi], N_proc, pflags);

    kalmanUpdate(Err[iP], Par[iP], msErr[hi], msPar[hi],
                 Err[iC], Par[iC], N_proc);

    if (Config::fit_val) MkFitter::CollectFitValidation(hi,N_proc,ev);
  }
  // XXXXX What's with chi2?
}

void MkFitter::CollectFitValidation(const int hi, const int N_proc, const Event * ev) const
{
  for (int n = 0; n < N_proc; ++n)
  {
    const float upt = 1.f/Par[iC](n,3,0);
    const FitVal tmpfitval
    {
      Par[iP].ConstAt(n,2,0),
      std::sqrt(Err[iP].ConstAt(n,2,2)),
      getPhi(Par[iP].ConstAt(n,0,0),Par[iP].ConstAt(n,1,0)),
      std::sqrt(getPhiErr2(Par[iP](n,0,0),Par[iP](n,1,0),Err[iP](n,0,0),Err[iP](n,1,1),Err[iP](n,0,1))),
      upt,
      std::sqrt(Err[iC](n,3,3))*upt*upt,
      Par[iC](n,4,0),
      std::sqrt(Err[iC](n,4,4)),
      getEta(Par[iC](n,5,0)),
      std::sqrt(Err[iC](n,5,5)/std::sin(Par[iC](n,5,0)))
    };
    
    ev->validation_.collectFitInfo(tmpfitval,Label(n,0,0),hi);
  }
}

void MkFitter::FitTracksSteered(const bool is_barrel[], const int N_proc, const Event * ev, const PropagationFlags pflags)
{
  // Fitting loop.

  dprintf("MkFitter::FitTracksSteered %d %d %d\n", is_barrel[0], is_barrel[1], is_barrel[2]);

  for (int hi = 0; hi < Nhits; ++hi)
  {
    // Note, charge is not passed (line propagation).
    // propagateLineToRMPlex(Err[iC], Par[iC], msErr[hi], msPar[hi],
    //                       Err[iP], Par[iP]);

    if (is_barrel[hi])
    {
      PropagateTracksToHitR(msPar[hi], N_proc, pflags);

      kalmanUpdate(Err[iP], Par[iP], msErr[hi], msPar[hi],
                   Err[iC], Par[iC], N_proc);
    }
    else
    {
      PropagateTracksToHitZ(msPar[hi], N_proc, pflags);

      kalmanUpdateEndcap(Err[iP], Par[iP], msErr[hi], msPar[hi],
                         Err[iC], Par[iC], N_proc);
    }

    if (Config::fit_val) MkFitter::CollectFitValidation(hi,N_proc,ev);
  }
  // XXXXX What's with chi2?
}

void MkFitter::OutputTracks(std::vector<Track>& tracks, int beg, int end, int iCP) const
{
  // Copies last track parameters (updated) into Track objects.
  // The tracks vector should be resized to allow direct copying.

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Err[iCP].CopyOut(itrack, tracks[i].errors_nc().Array());
    Par[iCP].CopyOut(itrack, tracks[i].parameters_nc().Array());
    
    tracks[i].setCharge(Chg(itrack, 0, 0));
    
    // XXXXX chi2 is not set (also not in SMatrix fit, it seems)
    tracks[i].setChi2(Chi2(itrack, 0, 0));
    tracks[i].setLabel(Label(itrack, 0, 0));
  }
}

void MkFitter::OutputFittedTracksAndHitIdx(std::vector<Track>& tracks, int beg, int end,
                                           bool outputProp) const
{
  // Copies last track parameters (updated) into Track objects and up to Nhits.
  // The tracks vector should be resized to allow direct copying.

  const int iO = outputProp ? iP : iC;

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Err[iO].CopyOut(itrack, tracks[i].errors_nc().Array());
    Par[iO].CopyOut(itrack, tracks[i].parameters_nc().Array());

    tracks[i].setCharge(Chg(itrack, 0, 0));
    tracks[i].setChi2(Chi2(itrack, 0, 0));
    tracks[i].setLabel(Label(itrack, 0, 0));

    tracks[i].resetHits();
    for (int hi = 0; hi < Nhits; ++hi)
    {
      tracks[i].addHitIdx(HoTArr[hi](itrack, 0, 0), 0.);
    }
  }
}

} // end namespace mkfit
