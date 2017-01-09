#include "MkFitter.h"
#include "CandCloner.h"

#include "PropagationMPlex.h"
#include "KalmanUtilsMPlex.h"
#include "ConformalUtilsMPlex.h"
#ifdef USE_CUDA
//#include "FitterCU.h"
#endif

//#define DEBUG
#include "Debug.h"

#include <sstream>

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
    printf("%5.2f  ", hipo(Par[idx].At(i, 3, 0), Par[idx].At(i, 4, 0)));
  }
}

int MkFitter::countValidHits(int itrack, int end_hit) const
{
  int result = 0;
  for (int hi = 0; hi < end_hit; ++hi)
    {
      if (HitsIdx[hi](itrack, 0, 0) >= 0) result++;
    }
  return result;
}

int MkFitter::countInvalidHits(int itrack, int end_hit) const
{
  int result = 0;
  for (int hi = 0; hi < end_hit; ++hi)
    {
      // XXXX MT: Should also count -2 hits as invalid?
      if (HitsIdx[hi](itrack, 0, 0) == -1) result++;
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

  int itrack;

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
// FIXME: uncomment when track building is ported to GPU.
#if 1
//#ifndef USE_CUDA
    for (int hi = 0; hi < Nhits; ++hi)
    {
      const int hidx = trk.getHitIdx(hi);
      const Hit &hit = layerHits[hi][hidx];

      HitsIdx[hi](itrack, 0, 0) = hidx;
      if (hidx<0) continue;

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
      const Hit &hit = layerHits[hi].m_hits[hidx];

      msErr[hi].CopyIn(itrack, hit.errArray());
      msPar[hi].CopyIn(itrack, hit.posArray());
      HitsIdx[hi](itrack, 0, 0) = hidx;
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

  const Track &trk = tracks[beg];
  const char *varr       = (char*) &trk;
  const int   off_error  = (char*) trk.errors().Array() - varr;
  const int   off_param  = (char*) trk.parameters().Array() - varr;

  int idx[NN]      __attribute__((aligned(64)));
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

    idx[itrack] = (char*) &trk - varr;

    Chg(itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();
  }

#ifdef MIC_INTRINSICS
  __m512i vi      = _mm512_load_epi32(idx);
  Err[iC].SlurpIn(varr + off_error, vi);
  Par[iC].SlurpIn(varr + off_param, vi);
#else
  Err[iC].SlurpIn(varr + off_error, idx);
  Par[iC].SlurpIn(varr + off_param, idx);
#endif
  
// CopyIn seems fast enough, but indirections are quite slow.
// For GPU computations, it has been moved in between kernels
// in an attempt to overlap CPU and GPU computations.
//#ifndef USE_CUDA
#if 1
  for (int hi = 0; hi < Nhits; ++hi)
  {
    const int   hidx      = tracks[beg].getHitIdx(hi);
    const Hit  &hit       = layerHits[hi][hidx];
    const char *varr      = (char*) &hit;
    const int   off_error = (char*) hit.errArray() - varr;
    const int   off_param = (char*) hit.posArray() - varr;

    for (int i = beg; i < end; ++i)
    {
      const int   hidx = tracks[i].getHitIdx(hi);
      const Hit  &hit  = layerHits[hi][hidx];
      itrack = i - beg;
      idx[itrack] = (char*) &hit - varr;
      HitsIdx[hi](itrack, 0, 0) = hidx;
    }

#ifdef MIC_INTRINSICS
    __m512i vi      = _mm512_load_epi32(idx);
    msErr[hi].SlurpIn(varr + off_error, vi);
    msPar[hi].SlurpIn(varr + off_param, vi);
#else
    msErr[hi].SlurpIn(varr + off_error, idx);
    msPar[hi].SlurpIn(varr + off_param, idx);
#endif
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

    Label(itrack, 0, 0) = trk.label();

    Err[iI].CopyIn(itrack, trk.errors().Array());
    Par[iI].CopyIn(itrack, trk.parameters().Array());

    Chg (itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();

    for (int hi = 0; hi < Nhits; ++hi)
    {
      // MPBL: It does not seem that these values are that dummies
      //       Not transfering them to the GPU reduces the number of
      //       nFoundHits in the printouts.
      HitsIdx[hi](itrack, 0, 0) = trk.getHitIdx(hi);//dummy value for now

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

      HitsIdx[hi](itrack, 0, 0) = trk.getHitIdx(hi);//dummy value for now

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
      const int hidx = trk.getHitIdx(hi);
      HitsIdx[hi](itrack, 0, 0) = hidx;
      if (hidx<0) continue;//fixme, check if this is harmless
      const Hit &hit = layerHits[hi][hidx];

      msErr[hi].CopyIn(itrack, hit.errArray());
      msPar[hi].CopyIn(itrack, hit.posArray());
    }
#endif
  }
}

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
#pragma simd
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

void MkFitter::FitTracks(const int N_proc, const Event * ev, const bool useParamBfield)
{
  // Fitting loop.

  for (int hi = 0; hi < Nhits; ++hi)
  {
    // Note, charge is not passed (line propagation).
    // propagateLineToRMPlex(Err[iC], Par[iC], msErr[hi], msPar[hi],
    //                       Err[iP], Par[iP]);

    propagateHelixToRMPlex(Err[iC], Par[iC], Chg, msPar[hi],
                           Err[iP], Par[iP], N_proc, useParamBfield);

    updateParametersMPlex(Err[iP], Par[iP], Chg, msErr[hi], msPar[hi],
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

void MkFitter::FitTracksTestEndcap(const int N_proc, const Event* ev, const bool useParamBfield)
{

  if (countValidHits(0)<8) return;
  if (Label.ConstAt(0, 0, 0)<0) return;

  //float ptin = 1./Par[iC].ConstAt(0, 0, 3);

  float chi2 = 0.;
  int hitcount = 0;
  for (int hi = 0; hi < Nhits; ++hi)
    {

      //if starting from seed, skip seed hits in fit (note there are only two hits in seeds since pxb1 is removed upfront, at least for now)
      if (Config::readCmsswSeeds && hi<2) continue;
      if (HitsIdx[hi].ConstAt(0, 0, 0)<0) continue;

      dprint("hit #" << hi << " hit  pos=" << msPar[hi].ConstAt(0, 0, 0) << ", " <<  msPar[hi].ConstAt(0, 0, 1) << ", " <<  msPar[hi].ConstAt(0, 0, 2));

      propagateHelixToZMPlex(Err[iC], Par[iC], Chg, msPar[hi],
			     Err[iP], Par[iP], N_proc, useParamBfield);

      dprint("hit #" << hi << " prop par=" << Par[iP].ConstAt(0, 0, 0) << ", " <<  Par[iP].ConstAt(0, 0, 1) << ", " <<  Par[iP].ConstAt(0, 0, 2)  << ", "
	     << Par[iP].ConstAt(0, 0, 3) << ", " <<  Par[iP].ConstAt(0, 0, 4) << ", " <<  Par[iP].ConstAt(0, 0, 5)
	     << " p=(" << cos(Par[iP].ConstAt(0, 0, 4))/Par[iP].ConstAt(0, 0, 3) << ", " << sin(Par[iP].ConstAt(0, 0, 4))/Par[iP].ConstAt(0, 0, 3)
	     << ", " << 1./(tan(Par[iP].ConstAt(0, 0, 5))*Par[iP].ConstAt(0, 0, 3)) << ")" );

      computeChi2EndcapMPlex(Err[iP], Par[iP], Chg, msErr[hi], msPar[hi],
			     Chi2,N_proc);

      dprint("hit chi2: " << Chi2.At(0, 0, 0));
      //if (Chi2.At(0, 0, 0)>30.) continue;
      chi2+=Chi2.At(0, 0, 0);

      updateParametersEndcapMPlex(Err[iP], Par[iP], Chg, msErr[hi], msPar[hi],
				  Err[iC], Par[iC], N_proc);

      dprint("hit #" << hi << " upda par=" << Par[iC].ConstAt(0, 0, 0) << ", " <<  Par[iC].ConstAt(0, 0, 1) << ", " <<  Par[iC].ConstAt(0, 0, 2)  << ", "
	     << Par[iC].ConstAt(0, 0, 3) << ", " <<  Par[iC].ConstAt(0, 0, 4) << ", " <<  Par[iC].ConstAt(0, 0, 5)
	     << " p=(" << cos(Par[iC].ConstAt(0, 0, 4))/Par[iC].ConstAt(0, 0, 3) << ", " << sin(Par[iC].ConstAt(0, 0, 4))/Par[iC].ConstAt(0, 0, 3)
	     << ", " << 1./(tan(Par[iC].ConstAt(0, 0, 5))*Par[iC].ConstAt(0, 0, 3)) << ")" );
      hitcount++;
    }

  // if ((1./Par[iC].ConstAt(0, 0, 3))<0. || chi2<0 || chi2>1000.) {
  //   printf("ANOMALY pt=%6.1f chi2=%6.1f seed pt=%6.1f q=%2i sim pt=%6.1f q=%2i flip=%2i \n",1./Par[iC].ConstAt(0, 0, 3),chi2,ptin,Chg.ConstAt(0, 0, 0),
  // 	   ev->simTracks_[Label.ConstAt(0, 0, 0)].pT(),ev->simTracks_[Label.ConstAt(0, 0, 0)].charge(),std::abs(Chg.ConstAt(0, 0, 0)-ev->simTracks_[Label.ConstAt(0, 0, 0)].charge())/2);
  // }

  // std::cout << "found track with pt: " << 1./Par[iC].ConstAt(0, 0, 3) << " chi2: " << chi2
  // 	    << " delta(pT)/sim_pT: " << ((1./Par[iC].ConstAt(0, 0, 3))-ev->simTracks_[Label.ConstAt(0, 0, 0)].pT())/ev->simTracks_[Label.ConstAt(0, 0, 0)].pT()
  // 	    << " fitted hits: " << hitcount<< std::endl;
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

    // XXXXX chi2 is not set (also not in SMatrix fit, it seems)

    tracks[i].resetHits();
    for (int hi = 0; hi < Nhits; ++hi)
    {
      tracks[i].addHitIdx(HitsIdx[hi](itrack, 0, 0),0.);
    }
  }
}

void MkFitter::PropagateTracksToR(float R, const int N_proc)
{
    propagateHelixToRMPlex(Err[iC], Par[iC], Chg, R,
                           Err[iP], Par[iP], N_proc);
}

void MkFitter::SelectHitIndices(const LayerOfHits &layer_of_hits, const int N_proc, bool dump)
{
  const int   iI = iP;
  const float nSigmaPhi = 3;
  const float nSigmaZ   = 3;

  // Vectorizing this makes it run slower!
  //#pragma ivdep
  //#pragma simd
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    XHitSize[itrack] = 0;

    float z, phi, dz, dphi;
    {
      const float x = Par[iI].ConstAt(itrack, 0, 0);
      const float y = Par[iI].ConstAt(itrack, 1, 0);

      const float r2 = x*x + y*y;

      z   = Par[iI].ConstAt(itrack, 2, 0);
      phi = getPhi(x, y);
      dz  = nSigmaZ * std::sqrt(Err[iI].ConstAt(itrack, 2, 2));

      const float dphidx = -y/r2, dphidy = x/r2;
      const float dphi2  = dphidx * dphidx * Err[iI].ConstAt(itrack, 0, 0) +
                           dphidy * dphidy * Err[iI].ConstAt(itrack, 1, 1) +
                       2 * dphidx * dphidy * Err[iI].ConstAt(itrack, 0, 1);

#ifdef HARD_CHECK
      assert(dphi2 >= 0);
#endif

      dphi = nSigmaPhi * std::sqrt(std::abs(dphi2));

      if (std::abs(dphi)<Config::minDPhi) dphi = Config::minDPhi;
      if (std::abs(dz)<Config::minDZ) dz = Config::minDZ;

      if (Config::useCMSGeom)
      {
        //now correct for bending and for layer thickness unsing linear approximation
        const float deltaR = Config::cmsDeltaRad; //fixme! using constant value, to be taken from layer properties
        const float r  = std::sqrt(r2);
#ifdef CCSCOORD
        //here alpha is the difference between posPhi and momPhi
        const float alpha = phi - Par[iP].ConstAt(itrack, 4, 0);
        float cosA, sinA;
        if (Config::useTrigApprox) {
          sincos4(alpha, sinA, cosA);
        } else {
          cosA = std::cos(alpha);
          sinA = std::sin(alpha);
        }
#else
        const float px = Par[iP].ConstAt(itrack, 3, 0);
        const float py = Par[iP].ConstAt(itrack, 4, 0);
        const float pt = std::sqrt(px*px + py*py);
        //here alpha is the difference between posPhi and momPhi
        const float cosA = ( x*px + y*py ) / (pt*r);
        const float sinA = ( y*px - x*py ) / (pt*r);
#endif
        //take abs so that we always inflate the window
        const float dist = std::abs(deltaR*sinA/cosA);
        dphi += dist / r;
      }
    }

    const LayerOfHits &L = layer_of_hits;

    if (std::abs(dz)   > Config::m_max_dz)   dz   = Config::m_max_dz;
    if (std::abs(dphi) > Config::m_max_dphi) dphi = Config::m_max_dphi;

    const int zb1 = L.GetZBinChecked(z - dz);
    const int zb2 = L.GetZBinChecked(z + dz) + 1;
    const int pb1 = L.GetPhiBin(phi - dphi);
    const int pb2 = L.GetPhiBin(phi + dphi) + 1;
    // MT: The extra phi bins give us ~1.5% more good tracks at expense of 10% runtime.
    // const int pb1 = L.GetPhiBin(phi - dphi) - 1;
    // const int pb2 = L.GetPhiBin(phi + dphi) + 2;

    if (dump)
      printf("LayerOfHits::SelectHitIndices %6.3f %6.3f %6.6f %7.5f %3d %3d %4d %4d\n",
             z, phi, dz, dphi, zb1, zb2, pb1, pb2);

    // MT: One could iterate in "spiral" order, to pick hits close to the center.
    // http://stackoverflow.com/questions/398299/looping-in-a-spiral
    // This would then work best with relatively small bin sizes.
    // Or, set them up so I can always take 3x3 array around the intersection.

    for (int zi = zb1; zi < zb2; ++zi)
    {
      for (int pi = pb1; pi < pb2; ++pi)
      {
        const int pb = pi & L.m_phi_mask;

        // MT: The following line is the biggest hog (4% total run time).
        // This comes from cache misses, I presume.
        // It might make sense to make first loop to extract bin indices
        // and issue prefetches at the same time.
        // Then enter vectorized loop to actually collect the hits in proper order.

        for (int hi = L.m_phi_bin_infos[zi][pb].first; hi < L.m_phi_bin_infos[zi][pb].second; ++hi)
        {
          // MT: Access into m_hit_zs and m_hit_phis is 1% run-time each.

#ifdef LOH_USE_PHI_Z_ARRAYS
          float ddz   = std::abs(z   - L.m_hit_zs[hi]);
          float ddphi = std::abs(phi - L.m_hit_phis[hi]);
          if (ddphi > Config::PI) ddphi = Config::TwoPI - ddphi;

          if (dump)
            printf("     SHI %3d %4d %4d %5d  %6.3f %6.3f %6.4f %7.5f   %s\n",
                   zi, pi, pb, hi,
                   L.m_hit_zs[hi], L.m_hit_phis[hi], ddz, ddphi,
                   (ddz < dz && ddphi < dphi) ? "PASS" : "FAIL");

          // MT: Commenting this check out gives full efficiency ...
          //     and means our error estimations are wrong!
          // Avi says we should have *minimal* search windows per layer.
          // Also ... if bins are sufficiently small, we do not need the extra
          // checks, see above.
          // if (ddz < dz && ddphi < dphi && XHitSize[itrack] < MPlexHitIdxMax)
#endif
          // MT: The following check also makes more sense with spiral traversal,
          // we'd be taking in closest hits first.
          if (XHitSize[itrack] < MPlexHitIdxMax)
          {
            XHitArr.At(itrack, XHitSize[itrack]++, 0) = hi;
          }
        }
      }
    }
  }
}

//==============================================================================
// AddBestHit()
//==============================================================================

//#define NO_PREFETCH
//#define NO_GATHER

void MkFitter::AddBestHit(const LayerOfHits &layer_of_hits, const int N_proc)
{
  float minChi2[NN];
  int   bestHit[NN];
  // MT: fill_n gave me crap on MIC, NN=8,16, doing in maxSize search below.
  // Must be a compiler issue.
  // std::fill_n(minChi2, NN, Config::chi2Cut);
  // std::fill_n(bestHit, NN, -1);

  const char *varr      = (char*) layer_of_hits.m_hits;

  const int   off_error = (char*) layer_of_hits.m_hits[0].errArray() - varr;
  const int   off_param = (char*) layer_of_hits.m_hits[0].posArray() - varr;

  int idx[NN]      __attribute__((aligned(64)));

  int maxSize = 0;

  // Determine maximum number of hits for tracks in the collection.
  // At the same time prefetch the first set of hits to L1 and the second one to L2.
  for (int it = 0; it < NN; ++it)
  {
    if (it < N_proc)
    {
      if (XHitSize[it] > 0)
      {
#ifndef NO_PREFETCH
        _mm_prefetch(varr + XHitArr.At(it, 0, 0) * sizeof(Hit), _MM_HINT_T0);
	if (XHitSize[it] > 1)
	{
	  _mm_prefetch(varr + XHitArr.At(it, 1, 0) * sizeof(Hit), _MM_HINT_T1);
        }
#endif
        maxSize = std::max(maxSize, XHitSize[it]);
      }
    }

    idx[it]     = 0;
    bestHit[it] = -1;
    minChi2[it] = Config::chi2Cut;
  }

// Has basically no effect, it seems.
//#pragma noprefetch
  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
    //fixme what if size is zero???

#pragma simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        idx[itrack] = XHitArr.At(itrack, hit_cnt, 0) * sizeof(Hit);
      }
    }
#if defined(MIC_INTRINSICS)
    __m512i vi = _mm512_load_epi32(idx);
#endif

#ifndef NO_PREFETCH
    // Prefetch to L2 the hits we'll process after two loops iterations.
    // Ideally this would be initiated before coming here, for whole bunch_of_hits.m_hits vector.
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt + 2 < XHitSize[itrack])
      {
        _mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+2, 0)*sizeof(Hit), _MM_HINT_T1);
      }
    }
#endif

#ifdef NO_GATHER

#pragma simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        const Hit &hit = layer_of_hits.m_hits[XHitArr.At(itrack, hit_cnt, 0)];
        msErr[Nhits].CopyIn(itrack, hit.errArray());
        msPar[Nhits].CopyIn(itrack, hit.posArray());
      }
    }
    
#else //NO_GATHER

#if defined(MIC_INTRINSICS)
    msErr[Nhits].SlurpIn(varr + off_error, vi);
    msPar[Nhits].SlurpIn(varr + off_param, vi);
#else
    msErr[Nhits].SlurpIn(varr + off_error, idx);
    msPar[Nhits].SlurpIn(varr + off_param, idx);
#endif
#endif //NO_GATHER

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    computeChi2MPlex(Err[iP], Par[iP], Chg, msErr[Nhits], msPar[Nhits], outChi2, N_proc);

#ifndef NO_PREFETCH
    // Prefetch to L1 the hits we'll process in the next loop iteration.
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt + 1 < XHitSize[itrack])
      {
        _mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+1, 0)*sizeof(Hit), _MM_HINT_T0);
      }
    }
#endif

    //update best hit in case chi2<minChi2
#pragma simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        const float chi2 = std::abs(outChi2[itrack]);//fixme negative chi2 sometimes...
        dprint("chi2=" << chi2 << " minChi2[itrack]=" << minChi2[itrack]);
        if (chi2 < minChi2[itrack])
        {
          minChi2[itrack] = chi2;
          bestHit[itrack] = XHitArr.At(itrack, hit_cnt, 0);
        }
      }
    }
  } // end loop over hits

  //copy in MkFitter the hit with lowest chi2
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    if (bestHit[itrack] >= 0)
    {
      _mm_prefetch( (const char*) & layer_of_hits.m_hits[ bestHit[itrack] ], _MM_HINT_T0);
    }
  }

#pragma simd
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    //fixme decide what to do in case no hit found
    if (bestHit[itrack] >= 0)
    {
      const Hit &hit  = layer_of_hits.m_hits[ bestHit[itrack] ];
      const float chi2 = minChi2[itrack];

      dprint("ADD BEST HIT FOR TRACK #" << itrack << std::endl
        << "prop x=" << Par[iP].ConstAt(itrack, 0, 0) << " y=" << Par[iP].ConstAt(itrack, 1, 0) << std::endl
        << "copy in hit #" << bestHit[itrack] << " x=" << hit.position()[0] << " y=" << hit.position()[1]);
	  
      msErr[Nhits].CopyIn(itrack, hit.errArray());
      msPar[Nhits].CopyIn(itrack, hit.posArray());
      Chi2(itrack, 0, 0) += chi2;
      HitsIdx[Nhits](itrack, 0, 0) = bestHit[itrack];
    }
    else
    {
      dprint("ADD FAKE HIT FOR TRACK #" << itrack);

      msErr[Nhits].SetDiagonal3x3(itrack, 666);
      msPar[Nhits](itrack,0,0) = Par[iP](itrack,0,0);
      msPar[Nhits](itrack,1,0) = Par[iP](itrack,1,0);
      msPar[Nhits](itrack,2,0) = Par[iP](itrack,2,0);
      HitsIdx[Nhits](itrack, 0, 0) = -1;

      // Don't update chi2
    }
  }

  //now update the track parameters with this hit (note that some calculations are already done when computing chi2... not sure it's worth caching them?)
  dprint("update parameters");
  updateParametersMPlex(Err[iP], Par[iP], Chg, msErr[Nhits], msPar[Nhits],
			Err[iC], Par[iC], N_proc);

  //std::cout << "Par[iP](0,0,0)=" << Par[iP](0,0,0) << " Par[iC](0,0,0)=" << Par[iC](0,0,0)<< std::endl;
}


void MkFitter::FindCandidates(const LayerOfHits &layer_of_hits,
                              std::vector<std::vector<Track> >& tmp_candidates,
                              const int offset, const int N_proc)
{
  const char *varr      = (char*) layer_of_hits.m_hits;

  const int   off_error = (char*) layer_of_hits.m_hits[0].errArray() - varr;
  const int   off_param = (char*) layer_of_hits.m_hits[0].posArray() - varr;

  int idx[NN]      __attribute__((aligned(64)));

  int maxSize = 0;

  // Determine maximum number of hits for tracks in the collection.
  // At the same time prefetch the first set of hits to L1 and the second one to L2.
  for (int it = 0; it < NN; ++it)
  {
    if (it < N_proc)
    {
      if (XHitSize[it] > 0)
      {
	_mm_prefetch(varr + XHitArr.At(it, 0, 0) * sizeof(Hit), _MM_HINT_T0);
	if (XHitSize[it] > 1)
	{
	  _mm_prefetch(varr + XHitArr.At(it, 1, 0) * sizeof(Hit), _MM_HINT_T1);
	}
	maxSize = std::max(maxSize, XHitSize[it]);
      }
    }

    idx[it] = 0;
    }

  // Has basically no effect, it seems.
  //#pragma noprefetch
  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
#pragma simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
	idx[itrack] = XHitArr.At(itrack, hit_cnt, 0) * sizeof(Hit);
      }
    }
#if defined(MIC_INTRINSICS)
    __m512i vi = _mm512_load_epi32(idx);
#endif
    
    // Prefetch to L2 the hits we'll (probably) process after two loops iterations.
    // Ideally this would be initiated before coming here, for whole bunch_of_hits.m_hits vector.
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt + 2 < XHitSize[itrack])
      {
	_mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+2, 0)*sizeof(Hit), _MM_HINT_T1);
      }
    }

#if defined(MIC_INTRINSICS)
    msErr[Nhits].SlurpIn(varr + off_error, vi);
    msPar[Nhits].SlurpIn(varr + off_param, vi);
#else
    msErr[Nhits].SlurpIn(varr + off_error, idx);
    msPar[Nhits].SlurpIn(varr + off_param, idx);
#endif

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    computeChi2MPlex(Err[iP], Par[iP], Chg, msErr[Nhits], msPar[Nhits], outChi2, N_proc);
    
    // Prefetch to L1 the hits we'll (probably) process in the next loop iteration.
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt + 1 < XHitSize[itrack])
      {
	_mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+1, 0)*sizeof(Hit), _MM_HINT_T0);
      }
    }
    
    //now update the track parameters with this hit (note that some calculations are already done when computing chi2, to be optimized)
    //this is not needed for candidates the hit is not added to, but it's vectorized so doing it serially below should take the same time
    //still it's a waste of time in case the hit is not added to any of the candidates, so check beforehand that at least one cand needs update
    bool oneCandPassCut = false;
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
	const float chi2 = std::abs(outChi2[itrack]);//fixme negative chi2 sometimes...
	dprint("chi2=" << chi2);
	if (chi2 < Config::chi2Cut)
	{
	  oneCandPassCut = true;
	  break;
	}
      }
    }
    
    if (oneCandPassCut)
    {
      updateParametersMPlex(Err[iP], Par[iP], Chg, msErr[Nhits], msPar[Nhits], Err[iC], Par[iC], N_proc);
      dprint("update parameters" << std::endl
	     << "propagated track parameters x=" << Par[iP].ConstAt(0, 0, 0) << " y=" << Par[iP].ConstAt(0, 1, 0) << std::endl
	     << "               hit position x=" << msPar[Nhits].ConstAt(0, 0, 0) << " y=" << msPar[Nhits].ConstAt(0, 1, 0) << std::endl
	     << "   updated track parameters x=" << Par[iC].ConstAt(0, 0, 0) << " y=" << Par[iC].ConstAt(0, 1, 0));
      
      //create candidate with hit in case chi2<Config::chi2Cut
      //fixme: please vectorize me... (not sure it's possible in this case)
      for (int itrack = 0; itrack < N_proc; ++itrack)
      {
	if (hit_cnt < XHitSize[itrack])
	{
	  const float chi2 = std::abs(outChi2[itrack]);//fixme negative chi2 sometimes...
	  dprint("chi2=" << chi2);
	  if (chi2 < Config::chi2Cut)
	  {
	    dprint("chi2 cut passed, creating new candidate");
	    //create a new candidate and fill the reccands_tmp vector
	    Track newcand;
	    newcand.resetHits();//probably not needed
	    newcand.setCharge(Chg(itrack, 0, 0));
	    newcand.setChi2(Chi2(itrack, 0, 0));
	    for (int hi = 0; hi < Nhits; ++hi)
	    {
	      newcand.addHitIdx(HitsIdx[hi](itrack, 0, 0),0.);//this should be ok since we already set the chi2 above
	    }
	    newcand.addHitIdx(XHitArr.At(itrack, hit_cnt, 0), chi2);
	    newcand.setLabel(Label(itrack, 0, 0));
	    //set the track state to the updated parameters
	    Err[iC].CopyOut(itrack, newcand.errors_nc().Array());
	    Par[iC].CopyOut(itrack, newcand.parameters_nc().Array());

	    dprint("updated track parameters x=" << newcand.parameters()[0] << " y=" << newcand.parameters()[1] << " z=" << newcand.parameters()[2] << " pt=" << 1./newcand.parameters()[3]);
	    
	    tmp_candidates[SeedIdx(itrack, 0, 0)-offset].push_back(newcand);
	  }
	}
      }
    }//end if (oneCandPassCut)
    
  }//end loop over hits
  
  //now add invalid hit
  //fixme: please vectorize me...
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    int hit_idx = countInvalidHits(itrack) < Config::maxHolesPerCand ? -1 : -2;
    Track newcand;
    newcand.resetHits();//probably not needed
    newcand.setCharge(Chg(itrack, 0, 0));
    newcand.setChi2(Chi2(itrack, 0, 0));
    for (int hi = 0; hi < Nhits; ++hi)
    {
      newcand.addHitIdx(HitsIdx[hi](itrack, 0, 0),0.);//this should be ok since we already set the chi2 above
    }
    newcand.addHitIdx(hit_idx, 0.);
    newcand.setLabel(Label(itrack, 0, 0));
    //set the track state to the propagated parameters
    Err[iP].CopyOut(itrack, newcand.errors_nc().Array());
    Par[iP].CopyOut(itrack, newcand.parameters_nc().Array());
    tmp_candidates[SeedIdx(itrack, 0, 0)-offset].push_back(newcand);
  }
}

void MkFitter::FindCandidatesEndcap(const LayerOfHits &layer_of_hits,
				    std::vector<std::vector<Track> >& tmp_candidates,
				    const int offset, const int N_proc)
{
  const char *varr      = (char*) layer_of_hits.m_hits;

  const int   off_error = (char*) layer_of_hits.m_hits[0].errArray() - varr;
  const int   off_param = (char*) layer_of_hits.m_hits[0].posArray() - varr;

  int idx[NN]      __attribute__((aligned(64)));

  int maxSize = 0;

  // Determine maximum number of hits for tracks in the collection.
  // At the same time prefetch the first set of hits to L1 and the second one to L2.
  for (int it = 0; it < NN; ++it)
  {
    if (it < N_proc)
    {
      if (XHitSize[it] > 0)
      {
	_mm_prefetch(varr + XHitArr.At(it, 0, 0) * sizeof(Hit), _MM_HINT_T0);
	if (XHitSize[it] > 1)
	{
	  _mm_prefetch(varr + XHitArr.At(it, 1, 0) * sizeof(Hit), _MM_HINT_T1);
	}
	maxSize = std::max(maxSize, XHitSize[it]);
      }
    }
    
    idx[it] = 0;
  }
  
  // Has basically no effect, it seems.
  //#pragma noprefetch
  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
#pragma simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
	idx[itrack] = XHitArr.At(itrack, hit_cnt, 0) * sizeof(Hit);
      }
    }
#if defined(MIC_INTRINSICS)
    __m512i vi = _mm512_load_epi32(idx);
#endif
      
    // Prefetch to L2 the hits we'll (probably) process after two loops iterations.
    // Ideally this would be initiated before coming here, for whole bunch_of_hits.m_hits vector.
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt + 2 < XHitSize[itrack])
      {
	_mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+2, 0)*sizeof(Hit), _MM_HINT_T1);
      }
    }
    
#if defined(MIC_INTRINSICS)    
    msErr[Nhits].SlurpIn(varr + off_error, vi);
    msPar[Nhits].SlurpIn(varr + off_param, vi);
#else
    msErr[Nhits].SlurpIn(varr + off_error, idx);
    msPar[Nhits].SlurpIn(varr + off_param, idx);
#endif

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    computeChi2EndcapMPlex(Err[iP], Par[iP], Chg, msErr[Nhits], msPar[Nhits], outChi2, N_proc);
    
    // Prefetch to L1 the hits we'll (probably) process in the next loop iteration.
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt + 1 < XHitSize[itrack])
      {
	_mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+1, 0)*sizeof(Hit), _MM_HINT_T0);
      }
    }
    
    //now update the track parameters with this hit (note that some calculations are already done when computing chi2, to be optimized)
    //this is not needed for candidates the hit is not added to, but it's vectorized so doing it serially below should take the same time
    //still it's a waste of time in case the hit is not added to any of the candidates, so check beforehand that at least one cand needs update
    bool oneCandPassCut = false;
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
	const float chi2 = std::abs(outChi2[itrack]);//fixme negative chi2 sometimes...
	dprint("chi2=" << chi2);
	if (chi2 < Config::chi2Cut)
	{
	  oneCandPassCut = true;
	  break;
	}
      }
    }
    
    if (oneCandPassCut)
    {
      updateParametersEndcapMPlex(Err[iP], Par[iP], Chg, msErr[Nhits], msPar[Nhits], Err[iC], Par[iC], N_proc);
      dprint("update parameters" << std::endl
	     << "propagated track parameters x=" << Par[iP].ConstAt(0, 0, 0) << " y=" << Par[iP].ConstAt(0, 1, 0) << " z=" << Par[iP].ConstAt(0, 2, 0) << " pt=" << 1./Par[iP].ConstAt(0, 3, 0) << std::endl
	     << "               hit position x=" << msPar[Nhits].ConstAt(0, 0, 0) << " y=" << msPar[Nhits].ConstAt(0, 1, 0) << std::endl
	     << "   updated track parameters x=" << Par[iC].ConstAt(0, 0, 0) << " y=" << Par[iC].ConstAt(0, 1, 0) << " z=" << Par[iC].ConstAt(0, 2, 0) << " pt=" << 1./Par[iC].ConstAt(0, 3, 0));
      
      //create candidate with hit in case chi2<Config::chi2Cut
      //fixme: please vectorize me... (not sure it's possible in this case)
      for (int itrack = 0; itrack < N_proc; ++itrack)
      {
	if (hit_cnt < XHitSize[itrack])
	{
	  const float chi2 = std::abs(outChi2[itrack]);//fixme negative chi2 sometimes...
	  dprint("chi2=" << chi2);
	  if (chi2 < Config::chi2Cut)
	  {
	    dprint("chi2 cut passed, creating new candidate");
	    //create a new candidate and fill the reccands_tmp vector
	    Track newcand;
	    newcand.resetHits();//probably not needed
	    newcand.setCharge(Chg(itrack, 0, 0));
	    newcand.setChi2(Chi2(itrack, 0, 0));
	    for (int hi = 0; hi < Nhits; ++hi)
	    {
	      newcand.addHitIdx(HitsIdx[hi](itrack, 0, 0),0.);//this should be ok since we already set the chi2 above
	    }
	    newcand.addHitIdx(XHitArr.At(itrack, hit_cnt, 0), chi2);
	    newcand.setLabel(Label(itrack, 0, 0));
	    //set the track state to the updated parameters
	    Err[iC].CopyOut(itrack, newcand.errors_nc().Array());
	    Par[iC].CopyOut(itrack, newcand.parameters_nc().Array());

	    dprint("updated track parameters x=" << newcand.parameters()[0] << " y=" << newcand.parameters()[1] << " z=" << newcand.parameters()[2] << " pt=" << 1./newcand.parameters()[3]);
	    
	    tmp_candidates[SeedIdx(itrack, 0, 0)-offset].push_back(newcand);
	  }
	}
      }
    }//end if (oneCandPassCut)
    
  }//end loop over hits
  
  //now add invalid hit
  //fixme: please vectorize me...
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    int hit_idx = countInvalidHits(itrack) < Config::maxHolesPerCand ? -1 : -2;
    
    bool withinBounds = true;
    float r2 = Par[iP](itrack,0,0)*Par[iP](itrack,0,0)+Par[iP](itrack,1,0)*Par[iP](itrack,1,0);
    if ( r2<(Config::cmsDiskMinRs[Nhits]*Config::cmsDiskMinRs[Nhits]) ||
	 r2>(Config::cmsDiskMaxRs[Nhits]*Config::cmsDiskMaxRs[Nhits]) ||
	 r2>(Config::cmsDiskMinRsHole[Nhits]*Config::cmsDiskMinRsHole[Nhits]) ||
	 r2<(Config::cmsDiskMaxRsHole[Nhits]*Config::cmsDiskMaxRsHole[Nhits]) ) withinBounds = false;
    //-3 means we did not expect any hit since we are out of bounds, so it does not count in countInvalidHits
    if (withinBounds==false) hit_idx = -3;

    Track newcand;
    newcand.resetHits();//probably not needed
    newcand.setCharge(Chg(itrack, 0, 0));
    newcand.setChi2(Chi2(itrack, 0, 0));
    for (int hi = 0; hi < Nhits; ++hi)
    {
      newcand.addHitIdx(HitsIdx[hi](itrack, 0, 0),0.);//this should be ok since we already set the chi2 above
    }
    newcand.addHitIdx(hit_idx, 0.);
    newcand.setLabel(Label(itrack, 0, 0));
      //set the track state to the propagated parameters
    Err[iP].CopyOut(itrack, newcand.errors_nc().Array());
    Par[iP].CopyOut(itrack, newcand.parameters_nc().Array());
    tmp_candidates[SeedIdx(itrack, 0, 0)-offset].push_back(newcand);
  }
}

void MkFitter::FindCandidatesMinimizeCopy(const LayerOfHits &layer_of_hits, CandCloner& cloner,
                                          const int offset, const int N_proc)
{
  const char *varr      = (char*) layer_of_hits.m_hits;

  const int   off_error = (char*) layer_of_hits.m_hits[0].errArray() - varr;
  const int   off_param = (char*) layer_of_hits.m_hits[0].posArray() - varr;

  int idx[NN]      __attribute__((aligned(64)));

  int maxSize = 0;

  // Determine maximum number of hits for tracks in the collection.
  // At the same time prefetch the first set of hits to L1 and the second one to L2.
#pragma simd
  for (int it = 0; it < NN; ++it)
  {
    if (it < N_proc)
    {
      if (XHitSize[it] > 0)
      {
        _mm_prefetch(varr + XHitArr.At(it, 0, 0) * sizeof(Hit), _MM_HINT_T0);
        if (XHitSize[it] > 1)
        {
          _mm_prefetch(varr + XHitArr.At(it, 1, 0) * sizeof(Hit), _MM_HINT_T1);
        }
        maxSize = std::max(maxSize, XHitSize[it]);
      }
    }

    idx[it] = 0;
  }
  // XXXX MT FIXME: use masks to filter out SlurpIns

// Has basically no effect, it seems.
//#pragma noprefetch
  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
#pragma simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        idx[itrack] = XHitArr.At(itrack, hit_cnt, 0) * sizeof(Hit);
      }
    }
#if defined(MIC_INTRINSICS)
    __m512i vi = _mm512_load_epi32(idx);
#endif

    // Prefetch to L2 the hits we'll (probably) process after two loops iterations.
    // Ideally this would be initiated before coming here, for whole bunch_of_hits.m_hits vector.
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt + 2 < XHitSize[itrack])
      {
        _mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+2, 0)*sizeof(Hit), _MM_HINT_T1);
      }
    }

#if defined(MIC_INTRINSICS)
    msErr[Nhits].SlurpIn(varr + off_error, vi);
    msPar[Nhits].SlurpIn(varr + off_param, vi);
#else
    msErr[Nhits].SlurpIn(varr + off_error, idx);
    msPar[Nhits].SlurpIn(varr + off_param, idx);
#endif

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    computeChi2MPlex(Err[iP], Par[iP], Chg, msErr[Nhits], msPar[Nhits], outChi2, N_proc);

    // Prefetch to L1 the hits we'll (probably) process in the next loop iteration.
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt + 1 < XHitSize[itrack])
      {
        _mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+1, 0)*sizeof(Hit), _MM_HINT_T0);
      }
    }

#pragma simd // DOES NOT VECTORIZE AS IT IS NOW
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      // make sure the hit was in the compatiblity window for the candidate

      if (hit_cnt < XHitSize[itrack])
      {
        float chi2 = fabs(outChi2[itrack]);//fixme negative chi2 sometimes...
#ifdef DEBUG
        std::cout << "chi2=" << chi2 << " for trkIdx=" << itrack << std::endl;
#endif
        if (chi2 < Config::chi2Cut)
        {
          IdxChi2List tmpList;
          tmpList.trkIdx = CandIdx(itrack, 0, 0);
          tmpList.hitIdx = XHitArr.At(itrack, hit_cnt, 0);
          tmpList.nhits  = countValidHits(itrack) + 1;
          tmpList.chi2   = Chi2(itrack, 0, 0) + chi2;
          cloner.add_cand(SeedIdx(itrack, 0, 0) - offset, tmpList);
          // hitsToAdd[SeedIdx(itrack, 0, 0)-offset].push_back(tmpList);
#ifdef DEBUG
          std::cout << "adding hit with hit_cnt=" << hit_cnt << " for trkIdx=" << tmpList.trkIdx << " orig Seed=" << Label(itrack, 0, 0) << std::endl;
#endif
        }
      }
    }

  }//end loop over hits

  //now add invalid hit
  //fixme: please vectorize me...
  for (int itrack = 0; itrack < N_proc; ++itrack)
    {
#ifdef DEBUG
      std::cout << "countInvalidHits(" << itrack << ")=" << countInvalidHits(itrack) << std::endl;
#endif

      int hit_idx = countInvalidHits(itrack) < Config::maxHolesPerCand ? -1 : -2;

      IdxChi2List tmpList;
      tmpList.trkIdx = CandIdx(itrack, 0, 0);
      tmpList.hitIdx = hit_idx;
      tmpList.nhits  = countValidHits(itrack);
      tmpList.chi2   = Chi2(itrack, 0, 0);
      cloner.add_cand(SeedIdx(itrack, 0, 0) - offset, tmpList);
      // hitsToAdd[SeedIdx(itrack, 0, 0)-offset].push_back(tmpList);
#ifdef DEBUG
      std::cout << "adding invalid hit" << std::endl;
#endif
    }
}



void MkFitter::InputTracksAndHitIdx(const std::vector<std::vector<Track> >& tracks,
                                    const std::vector<std::pair<int,MkFitter::IdxChi2List> >& idxs,
                                    int beg, int end, bool inputProp)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  //fixme: why do we need both i and itrack in the loops below?

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {

    const Track &trk = tracks[idxs[i].first][idxs[i].second.trkIdx];

    Label(itrack, 0, 0) = trk.label();
    SeedIdx(itrack, 0, 0) = idxs[i].first;
    CandIdx(itrack, 0, 0) = idxs[i].second.trkIdx;

    if (inputProp) 
    {
      Err[iP].CopyIn(itrack, trk.errors().Array());
      Par[iP].CopyIn(itrack, trk.parameters().Array());
    } 
    else 
    {      
      Err[iC].CopyIn(itrack, trk.errors().Array());
      Par[iC].CopyIn(itrack, trk.parameters().Array());
    } 

    Chg(itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();

    for (int hi = 0; hi < Nhits; ++hi)
    {
      HitsIdx[hi](itrack, 0, 0) = trk.getHitIdx(hi);//dummy value for now
    }
  }
}


void MkFitter::UpdateWithLastHit(const LayerOfHits &layer_of_hits, int N_proc)
{
  for (int i = 0; i < N_proc; ++i)
  {
    int hit_idx = HitsIdx[Nhits - 1](i, 0, 0);

    if (hit_idx < 0) continue;

    Hit &hit = layer_of_hits.m_hits[hit_idx];

    msErr[Nhits - 1].CopyIn(i, hit.errArray());
    msPar[Nhits - 1].CopyIn(i, hit.posArray());
  }

  updateParametersMPlex(Err[iP], Par[iP], Chg, msErr[Nhits-1], msPar[Nhits-1], Err[iC], Par[iC], N_proc);

  //now that we have moved propagation at the end of the sequence we lost the handle of
  //using the propagated parameters instead of the updated for the missing hit case.
  //so we need to replace by hand the updated with the propagated
  //there may be a better way to restore this...

  for (int i = 0; i < N_proc; ++i)
  {
    if (HitsIdx[Nhits - 1](i, 0, 0) < 0)
    {
      float tmp[21];
      Err[iP].CopyOut(i, tmp);
      Err[iC].CopyIn(i, tmp);
      Par[iP].CopyOut(i, tmp);
      Par[iC].CopyIn(i, tmp);
    }
  }
}

void MkFitter::UpdateWithLastHitEndcap(const LayerOfHits &layer_of_hits, int N_proc)
{
  for (int i = 0; i < N_proc; ++i)
  {
    int hit_idx = HitsIdx[Nhits - 1](i, 0, 0);

    if (hit_idx < 0) continue;

    Hit &hit = layer_of_hits.m_hits[hit_idx];

    msErr[Nhits - 1].CopyIn(i, hit.errArray());
    msPar[Nhits - 1].CopyIn(i, hit.posArray());
  }

  updateParametersEndcapMPlex(Err[iP], Par[iP], Chg, msErr[Nhits-1], msPar[Nhits-1], Err[iC], Par[iC], N_proc);

  //now that we have moved propagation at the end of the sequence we lost the handle of
  //using the propagated parameters instead of the updated for the missing hit case.
  //so we need to replace by hand the updated with the propagated
  //there may be a better way to restore this...

  for (int i = 0; i < N_proc; ++i)
  {
    if (HitsIdx[Nhits - 1](i, 0, 0) < 0)
    {
      float tmp[21];
      Err[iP].CopyOut(i, tmp);
      Err[iC].CopyIn(i, tmp);
      Par[iP].CopyOut(i, tmp);
      Par[iC].CopyIn(i, tmp);
    }
  }
}

void MkFitter::CopyOutParErr(std::vector<std::vector<Track> >& seed_cand_vec,
                             int N_proc, bool outputProp) const
{
  const int iO = outputProp ? iP : iC;

  for (int i = 0; i < N_proc; ++i)
  {
    //create a new candidate and fill the cands_for_next_lay vector
    Track &cand = seed_cand_vec[SeedIdx(i, 0, 0)][CandIdx(i, 0, 0)];

    //set the track state to the updated parameters
    Err[iO].CopyOut(i, cand.errors_nc().Array());
    Par[iO].CopyOut(i, cand.parameters_nc().Array());

    dprint((outputProp?"propagated":"updated") << " track parameters x=" << cand.parameters()[0]
              << " y=" << cand.parameters()[1]
              << " z=" << cand.parameters()[2]
              << " pt=" << 1./cand.parameters()[3]
              << " posEta=" << cand.posEta()
              << " etaBin=" << getEtaBin(cand.posEta()));
  }
}



//-----------------------------------------------------------------------------------------//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//-----------------------------------------------------------------------------------------//
// Endcap section: duplicate functions for now, cleanup and merge once strategy sorted out //
//-----------------------------------------------------------------------------------------//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//-----------------------------------------------------------------------------------------//

void MkFitter::PropagateTracksToZ(float Z, const int N_proc)
{
    propagateHelixToZMPlex(Err[iC], Par[iC], Chg, Z,
                           Err[iP], Par[iP], N_proc);
}

void MkFitter::SelectHitIndicesEndcap(const LayerOfHits &layer_of_hits, const int N_proc, bool dump)
{
  const int   iI = iP;
  const float nSigmaPhi = 3;
  const float nSigmaR   = 3;

  //dump = true;

  // Vectorizing this makes it run slower!
  //#pragma ivdep
  //#pragma simd
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    XHitSize[itrack] = 0;

    float r, phi, dr, dphi;
    {
      const float x = Par[iI].ConstAt(itrack, 0, 0);
      const float y = Par[iI].ConstAt(itrack, 1, 0);

      const float r2 = x*x + y*y;
      r  = std::sqrt(r2);

      phi = getPhi(x, y);
      dr  = nSigmaR*(x*x*Err[iI].ConstAt(itrack, 0, 0) + y*y*Err[iI].ConstAt(itrack, 1, 1) + 2*x*y*Err[iI].ConstAt(itrack, 0, 1))/r2;

      const float dphidx = -y/r2, dphidy = x/r2;
      const float dphi2  = dphidx * dphidx * Err[iI].ConstAt(itrack, 0, 0) +
                           dphidy * dphidy * Err[iI].ConstAt(itrack, 1, 1) +
                       2 * dphidx * dphidy * Err[iI].ConstAt(itrack, 0, 1);

#ifdef HARD_CHECK
      assert(dphi2 >= 0);
#endif

      dphi = nSigmaPhi * std::sqrt(std::abs(dphi2));

      if (std::abs(dphi)<Config::minDPhi) dphi = Config::minDPhi;
//       if (std::abs(dz)<Config::minDZ) dz = Config::minDZ;

      if (Config::useCMSGeom)
      {
#ifdef CCSCOORD
        //now correct for bending and for layer thickness unsing linear approximation
        const float deltaZ = 5; //fixme! using constant value, to be taken from layer properties
        float cosT = std::cos(Par[iI].ConstAt(itrack, 5, 0));
        float sinT = std::sin(Par[iI].ConstAt(itrack, 5, 0));
        //here alpha is the helix angular path corresponding to deltaZ
        const float k = Chg.ConstAt(itrack, 0, 0) * 100.f / (-Config::sol*Config::Bfield);
        const float alpha  = deltaZ*sinT*Par[iI].ConstAt(itrack, 3, 0)/(cosT*k);
        dphi += std::abs(alpha);
#else
        assert(0);
#endif
      }
    }

    const LayerOfHits &L = layer_of_hits;

    // if (std::abs(dz)   > Config::m_max_dz)   dz   = Config::m_max_dz;
    if (std::abs(dphi) > Config::m_max_dphi) dphi = Config::m_max_dphi;

    const int rb1 = L.GetRBinChecked(r - dr);
    const int rb2 = L.GetRBinChecked(r + dr) + 1;
    const int pb1 = L.GetPhiBin(phi - dphi);
    const int pb2 = L.GetPhiBin(phi + dphi) + 1;
    // MT: The extra phi bins give us ~1.5% more good tracks at expense of 10% runtime.
    // const int pb1 = L.GetPhiBin(phi - dphi) - 1;
    // const int pb2 = L.GetPhiBin(phi + dphi) + 2;

    if (dump)
      printf("LayerOfHits::SelectHitIndices %6.3f %6.3f %6.4f %7.5f %3d %3d %4d %4d\n",
             r, phi, dr, dphi, rb1, rb2, pb1, pb2);

    // MT: One could iterate in "spiral" order, to pick hits close to the center.
    // http://stackoverflow.com/questions/398299/looping-in-a-spiral
    // This would then work best with relatively small bin sizes.
    // Or, set them up so I can always take 3x3 array around the intersection.

    for (int ri = rb1; ri < rb2; ++ri)
    {
      for (int pi = pb1; pi < pb2; ++pi)
      {
        const int pb = pi & L.m_phi_mask;

        // MT: The following line is the biggest hog (4% total run time).
        // This comes from cache misses, I presume.
        // It might make sense to make first loop to extract bin indices
        // and issue prefetches at the same time.
        // Then enter vectorized loop to actually collect the hits in proper order.

        for (int hi = L.m_phi_bin_infos[ri][pb].first; hi < L.m_phi_bin_infos[ri][pb].second; ++hi)
        {
          // MT: Access into m_hit_zs and m_hit_phis is 1% run-time each.

#ifdef LOH_USE_PHI_Z_ARRAYS
          float ddz   = std::abs(z   - L.m_hit_zs[hi]);
          float ddphi = std::abs(phi - L.m_hit_phis[hi]);
          if (ddphi > Config::PI) ddphi = Config::TwoPI - ddphi;

          if (dump)
            printf("     SHI %3d %4d %4d %5d  %6.3f %6.3f %6.4f %7.5f   %s\n",
                   ri, pi, pb, hi,
                   L.m_hit_zs[hi], L.m_hit_phis[hi], ddz, ddphi,
                   (ddz < dz && ddphi < dphi) ? "PASS" : "FAIL");

          // MT: Commenting this check out gives full efficiency ...
          //     and means our error estimations are wrong!
          // Avi says we should have *minimal* search windows per layer.
          // Also ... if bins are sufficiently small, we do not need the extra
          // checks, see above.
          // if (ddz < dz && ddphi < dphi && XHitSize[itrack] < MPlexHitIdxMax)
#endif
          // MT: The following check also makes more sense with spiral traversal,
          // we'd be taking in closest hits first.
          if (XHitSize[itrack] < MPlexHitIdxMax)
          {
            XHitArr.At(itrack, XHitSize[itrack]++, 0) = hi;
          }
        }
      }
    }
  }
}

void MkFitter::AddBestHitEndcap(const LayerOfHits &layer_of_hits, const int N_proc)
{
  float minChi2[NN];
  int   bestHit[NN];
  // MT: fill_n gave me crap on MIC, NN=8,16, doing in maxSize search below.
  // Must be a compiler issue.
  // std::fill_n(minChi2, NN, Config::chi2Cut);
  // std::fill_n(bestHit, NN, -1);

  const char *varr      = (char*) layer_of_hits.m_hits;

  const int   off_error = (char*) layer_of_hits.m_hits[0].errArray() - varr;
  const int   off_param = (char*) layer_of_hits.m_hits[0].posArray() - varr;

  int idx[NN]      __attribute__((aligned(64)));

  int maxSize = 0;

  // Determine maximum number of hits for tracks in the collection.
  // At the same time prefetch the first set of hits to L1 and the second one to L2.
  for (int it = 0; it < NN; ++it)
  {
    if (it < N_proc)
    {
      if (XHitSize[it] > 0)
      {
#ifndef NO_PREFETCH
        _mm_prefetch(varr + XHitArr.At(it, 0, 0) * sizeof(Hit), _MM_HINT_T0);
	if (XHitSize[it] > 1)
	{
	  _mm_prefetch(varr + XHitArr.At(it, 1, 0) * sizeof(Hit), _MM_HINT_T1);
        }
#endif
        maxSize = std::max(maxSize, XHitSize[it]);
      }
    }

    idx[it]     = 0;
    bestHit[it] = -1;
    minChi2[it] = Config::chi2Cut;
  }

// Has basically no effect, it seems.
//#pragma noprefetch
  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
    //fixme what if size is zero???

#pragma simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        idx[itrack] = XHitArr.At(itrack, hit_cnt, 0) * sizeof(Hit);
      }
    }
#if defined(MIC_INTRINSICS)
    __m512i vi = _mm512_load_epi32(idx);
#endif

#ifndef NO_PREFETCH
    // Prefetch to L2 the hits we'll process after two loops iterations.
    // Ideally this would be initiated before coming here, for whole bunch_of_hits.m_hits vector.
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt + 2 < XHitSize[itrack])
      {
        _mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+2, 0)*sizeof(Hit), _MM_HINT_T1);
      }
    }
#endif

#ifdef NO_GATHER

#pragma simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        const Hit &hit = layer_of_hits.m_hits[XHitArr.At(itrack, hit_cnt, 0)];
        msErr[Nhits].CopyIn(itrack, hit.errArray());
        msPar[Nhits].CopyIn(itrack, hit.posArray());
      }
    }

#else //NO_GATHER

#if defined(MIC_INTRINSICS)
    msErr[Nhits].SlurpIn(varr + off_error, vi);
    msPar[Nhits].SlurpIn(varr + off_param, vi);
#else
    msErr[Nhits].SlurpIn(varr + off_error, idx);
    msPar[Nhits].SlurpIn(varr + off_param, idx);
#endif
#endif //NO_GATHER

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    computeChi2EndcapMPlex(Err[iP], Par[iP], Chg, msErr[Nhits], msPar[Nhits], outChi2, N_proc);

#ifndef NO_PREFETCH
    // Prefetch to L1 the hits we'll process in the next loop iteration.
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt + 1 < XHitSize[itrack])
      {
        _mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+1, 0)*sizeof(Hit), _MM_HINT_T0);
      }
    }
#endif

    //update best hit in case chi2<minChi2
#pragma simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        const float chi2 = std::abs(outChi2[itrack]);//fixme negative chi2 sometimes...
        dprint("chi2=" << chi2 << " minChi2[itrack]=" << minChi2[itrack]);
        if (chi2 < minChi2[itrack])
        {
          minChi2[itrack] = chi2;
          bestHit[itrack] = XHitArr.At(itrack, hit_cnt, 0);
        }
      }
    }
  } // end loop over hits

  //copy in MkFitter the hit with lowest chi2
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    if (bestHit[itrack] >= 0)
    {
      _mm_prefetch( (const char*) & layer_of_hits.m_hits[ bestHit[itrack] ], _MM_HINT_T0);
    }
  }

#pragma simd
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    //fixme decide what to do in case no hit found
    if (bestHit[itrack] >= 0)
    {
      const Hit &hit  = layer_of_hits.m_hits[ bestHit[itrack] ];
      const float chi2 = minChi2[itrack];

      dprint("ADD BEST HIT FOR TRACK #" << itrack << std::endl
        << "prop x=" << Par[iP].ConstAt(itrack, 0, 0) << " y=" << Par[iP].ConstAt(itrack, 1, 0) << std::endl
        << "copy in hit #" << bestHit[itrack] << " x=" << hit.position()[0] << " y=" << hit.position()[1]);

      msErr[Nhits].CopyIn(itrack, hit.errArray());
      msPar[Nhits].CopyIn(itrack, hit.posArray());
      Chi2(itrack, 0, 0) += chi2;
      HitsIdx[Nhits](itrack, 0, 0) = bestHit[itrack];
    }
    else
    {

      bool withinBounds = true;
      float r2 = Par[iP](itrack,0,0)*Par[iP](itrack,0,0)+Par[iP](itrack,1,0)*Par[iP](itrack,1,0);
      if ( r2<(Config::cmsDiskMinRs[Nhits]*Config::cmsDiskMinRs[Nhits]) ||
	   r2>(Config::cmsDiskMaxRs[Nhits]*Config::cmsDiskMaxRs[Nhits]) ||
	   r2>(Config::cmsDiskMinRsHole[Nhits]*Config::cmsDiskMinRsHole[Nhits]) ||
	   r2<(Config::cmsDiskMaxRsHole[Nhits]*Config::cmsDiskMaxRsHole[Nhits]) ) withinBounds = false;

#ifdef DEBUG
      r2= std::sqrt(r2);
#endif
      dprint("ADD FAKE HIT FOR TRACK #" << itrack << " withinBounds=" << withinBounds << " r=" << r2
	     << " minR=" << Config::cmsDiskMinRs[Nhits] << " maxR=" << Config::cmsDiskMaxRs[Nhits]
	     << " minRHole=" << Config::cmsDiskMinRsHole[Nhits] << " maxRHole=" << Config::cmsDiskMaxRsHole[Nhits]);

      msErr[Nhits].SetDiagonal3x3(itrack, 666);
      msPar[Nhits](itrack,0,0) = Par[iP](itrack,0,0);
      msPar[Nhits](itrack,1,0) = Par[iP](itrack,1,0);
      msPar[Nhits](itrack,2,0) = Par[iP](itrack,2,0);
      //-3 means we did not expect any hit since we are out of bounds, so it does not count in countInvalidHits
      HitsIdx[Nhits](itrack, 0, 0) = withinBounds ? -1 : -3;

      // Don't update chi2
    }
  }

  //now update the track parameters with this hit (note that some calculations are already done when computing chi2... not sure it's worth caching them?)
  dprint("update parameters");
  updateParametersEndcapMPlex(Err[iP], Par[iP], Chg, msErr[Nhits], msPar[Nhits],
			      Err[iC], Par[iC], N_proc);

  //std::cout << "Par[iP](0,0,0)=" << Par[iP](0,0,0) << " Par[iC](0,0,0)=" << Par[iC](0,0,0)<< std::endl;
}

void MkFitter::FindCandidatesMinimizeCopyEndcap(const LayerOfHits &layer_of_hits, CandCloner& cloner,
                                                const int offset, const int N_proc)
{
  const char *varr      = (char*) layer_of_hits.m_hits;

  const int   off_error = (char*) layer_of_hits.m_hits[0].errArray() - varr;
  const int   off_param = (char*) layer_of_hits.m_hits[0].posArray() - varr;

  int idx[NN]      __attribute__((aligned(64)));

  int maxSize = 0;

  // Determine maximum number of hits for tracks in the collection.
  // At the same time prefetch the first set of hits to L1 and the second one to L2.
#pragma simd
  for (int it = 0; it < NN; ++it)
  {
    if (it < N_proc)
    {
      if (XHitSize[it] > 0)
      {
        _mm_prefetch(varr + XHitArr.At(it, 0, 0) * sizeof(Hit), _MM_HINT_T0);
        if (XHitSize[it] > 1)
        {
          _mm_prefetch(varr + XHitArr.At(it, 1, 0) * sizeof(Hit), _MM_HINT_T1);
        }
        maxSize = std::max(maxSize, XHitSize[it]);
      }
    }

    idx[it] = 0;
  }
  // XXXX MT FIXME: use masks to filter out SlurpIns

// Has basically no effect, it seems.
//#pragma noprefetch
  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt)
  {
#pragma simd
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        idx[itrack] = XHitArr.At(itrack, hit_cnt, 0) * sizeof(Hit);
      }
    }
#if defined(MIC_INTRINSICS)
    __m512i vi = _mm512_load_epi32(idx);
#endif

    // Prefetch to L2 the hits we'll (probably) process after two loops iterations.
    // Ideally this would be initiated before coming here, for whole bunch_of_hits.m_hits vector.
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt + 2 < XHitSize[itrack])
      {
        _mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+2, 0)*sizeof(Hit), _MM_HINT_T1);
      }
    }

#if defined(MIC_INTRINSICS)
    msErr[Nhits].SlurpIn(varr + off_error, vi);
    msPar[Nhits].SlurpIn(varr + off_param, vi);
#else
    msErr[Nhits].SlurpIn(varr + off_error, idx);
    msPar[Nhits].SlurpIn(varr + off_param, idx);
#endif

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    computeChi2EndcapMPlex(Err[iP], Par[iP], Chg, msErr[Nhits], msPar[Nhits], outChi2, N_proc);

    // Prefetch to L1 the hits we'll (probably) process in the next loop iteration.
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      if (hit_cnt + 1 < XHitSize[itrack])
      {
        _mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+1, 0)*sizeof(Hit), _MM_HINT_T0);
      }
    }

#pragma simd // DOES NOT VECTORIZE AS IT IS NOW
    for (int itrack = 0; itrack < N_proc; ++itrack)
    {
      // make sure the hit was in the compatiblity window for the candidate

      if (hit_cnt < XHitSize[itrack])
      {
        float chi2 = fabs(outChi2[itrack]);//fixme negative chi2 sometimes...
#ifdef DEBUG
        std::cout << "chi2=" << chi2 << " for trkIdx=" << itrack << std::endl;
#endif
        if (chi2 < Config::chi2Cut)
        {
          IdxChi2List tmpList;
          tmpList.trkIdx = CandIdx(itrack, 0, 0);
          tmpList.hitIdx = XHitArr.At(itrack, hit_cnt, 0);
          tmpList.nhits  = countValidHits(itrack) + 1;
          tmpList.chi2   = Chi2(itrack, 0, 0) + chi2;
          cloner.add_cand(SeedIdx(itrack, 0, 0) - offset, tmpList);
          // hitsToAdd[SeedIdx(itrack, 0, 0)-offset].push_back(tmpList);
#ifdef DEBUG
          std::cout << "adding hit with hit_cnt=" << hit_cnt << " for trkIdx=" << tmpList.trkIdx << " orig Seed=" << Label(itrack, 0, 0) << std::endl;
#endif
        }
      }
    }

  }//end loop over hits

  //now add invalid hit
  //fixme: please vectorize me...
  for (int itrack = 0; itrack < N_proc; ++itrack)
    {
#ifdef DEBUG
      std::cout << "countInvalidHits(" << itrack << ")=" << countInvalidHits(itrack) << std::endl;
#endif

      int hit_idx = countInvalidHits(itrack) < Config::maxHolesPerCand ? -1 : -2;

      bool withinBounds = true;
      float r2 = Par[iP](itrack,0,0)*Par[iP](itrack,0,0)+Par[iP](itrack,1,0)*Par[iP](itrack,1,0);
      if ( r2<(Config::cmsDiskMinRs[Nhits]*Config::cmsDiskMinRs[Nhits]) ||
           r2>(Config::cmsDiskMaxRs[Nhits]*Config::cmsDiskMaxRs[Nhits]) ||
           r2>(Config::cmsDiskMinRsHole[Nhits]*Config::cmsDiskMinRsHole[Nhits]) ||
           r2<(Config::cmsDiskMaxRsHole[Nhits]*Config::cmsDiskMaxRsHole[Nhits]) ) withinBounds = false;
      //-3 means we did not expect any hit since we are out of bounds, so it does not count in countInvalidHits
      if (withinBounds==false) hit_idx = -3;

      IdxChi2List tmpList;
      tmpList.trkIdx = CandIdx(itrack, 0, 0);
      tmpList.hitIdx = hit_idx;
      tmpList.nhits  = countValidHits(itrack);
      tmpList.chi2   = Chi2(itrack, 0, 0);
      cloner.add_cand(SeedIdx(itrack, 0, 0) - offset, tmpList);
      // hitsToAdd[SeedIdx(itrack, 0, 0)-offset].push_back(tmpList);
#ifdef DEBUG
      std::cout << "adding invalid hit" << std::endl;
#endif
    }
}
