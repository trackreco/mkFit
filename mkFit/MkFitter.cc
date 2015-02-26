#include "MkFitter.h"

#include "PropagationMPlex.h"
#include "KalmanUtilsMPlex.h"

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

//==============================================================================

void MkFitter::InputTracksAndHits(std::vector<Track>& tracks, int beg, int end)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {

    Track &trk = tracks[i];

    Err[iC].CopyIn(itrack, trk.errors().Array());
    Par[iC].CopyIn(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();

    for (int hi = 0; hi < Nhits; ++hi)
    {

      Hit &hit = trk.hitsVector()[hi];

      msErr[hi].CopyIn(itrack, hit.error().Array());
      msPar[hi].CopyIn(itrack, hit.parameters().Array());
    }
  }
}

void MkFitter::InputTracksOnly(std::vector<Track>& tracks, int beg, int end)
{
  // Assign track parameters to initial state, do NOT copy hit values.
  // Used for benchmarking the fitting with less "copy-in" load.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Track &trk = tracks[i];

    Err[iC].CopyIn(itrack, trk.errors().Array());
    Par[iC].CopyIn(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();
  }
}

void MkFitter::InputHitsOnly(std::vector<Hit>& hits, int beg, int end)
{
  // Push hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Hit &hit = hits[itrack];

    msErr[Nhits].CopyIn(itrack, hit.error().Array());
    msPar[Nhits].CopyIn(itrack, hit.parameters().Array());
  }
  Nhits++;
}

void MkFitter::FitTracks()
{
  // Fitting loop.

  for (int hi = 0; hi < Nhits; ++hi)
  {
    // Note, charge is not passed (line propagation).
    // propagateLineToRMPlex(Err[iC], Par[iC], msErr[hi], msPar[hi],
    //                       Err[iP], Par[iP]);

    propagateHelixToRMPlex(Err[iC], Par[iC], Chg, msPar[hi],
                           Err[iP], Par[iP]);

    updateParametersMPlex(Err[iP], Par[iP], msErr[hi], msPar[hi],
                          Err[iC], Par[iC]);
  }

  // XXXXX What's with chi2?
}

void MkFitter::OutputTracks(std::vector<Track>& tracks, int beg, int end, int iCP)
{
  // Copies last track parameters (updated) into Track objects.
  // The tracks vector should be resized to allow direct copying.

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Err[iCP].CopyOut(itrack, tracks[i].errors().Array());
    Par[iCP].CopyOut(itrack, tracks[i].parameters().Array());

    tracks[i].charge() = Chg(itrack, 0, 0);

    // XXXXX chi2 is not set (also not in SMatrix fit, it seems)
    tracks[i].setChi2(Chi2(itrack, 0, 0));
  }
}

void MkFitter::OutputFittedTracksAndHits(std::vector<Track>& tracks, int beg, int end)
{
  // Copies last track parameters (updated) into Track objects and up to Nhits.
  // The tracks vector should be resized to allow direct copying.

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Err[iC].CopyOut(itrack, tracks[i].errors().Array());
    Par[iC].CopyOut(itrack, tracks[i].parameters().Array());

    tracks[i].charge() = Chg(itrack, 0, 0);
    tracks[i].setChi2(Chi2(itrack, 0, 0));

    // XXXXX chi2 is not set (also not in SMatrix fit, it seems)

    tracks[i].resetHits();
    for (int hi = 0; hi < Nhits; ++hi)
    {

      Hit hit;

      msErr[hi].CopyOut(itrack, hit.error().Array());
      msPar[hi].CopyOut(itrack, hit.parameters().Array());

      tracks[i].addHit(hit,0.);

    }

  }
}

void MkFitter::PropagateTracksToR(float R)
{

    propagateHelixToRMPlex(Err[iC], Par[iC], Chg, R,
                           Err[iP], Par[iP]);

}

//fixme: do it properly with phi segmentation
void MkFitter::AddBestHit(std::vector<Hit>& lay_hits, int beg, int end)
{

  //fixme solve ambiguity NN vs beg-end
  float minChi2[NN];
  std::fill_n(minChi2, NN, 9999.);
  int bestHit[NN];
  std::fill_n(bestHit, NN, -1);

  //outer loop over hits, so that tracks can be vectorized
  int ih = 0;
  for ( ; ih<lay_hits.size(); ++ih)
  {

    std::cout << "consider hit #" << ih << " of " << lay_hits.size() << std::endl;
    Hit &hit = lay_hits[ih];
    std::cout << "hit x=" << hit.position()[0] << " y=" << hit.position()[1] << std::endl;      

    //create a dummy matriplex with N times the same hit
    //fixme this is just because I'm lazy and do not want to rewrite the chi2 calculation for simple vectors mixed with matriplex
    MPlexHS msErr_oneHit;
    MPlexHV msPar_oneHit;
    int itrack = 0;
    //fixme: please vectorize me...
    for (int i = beg; i < end; ++i, ++itrack)
    {
      msErr_oneHit.CopyIn(itrack, hit.error().Array());
      msPar_oneHit.CopyIn(itrack, hit.parameters().Array());
    }

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    computeChi2MPlex(Err[iP], Par[iP],msErr_oneHit, msPar_oneHit, outChi2);

    //update best hit in case chi2<minChi2
    itrack = 0;
    //fixme: please vectorize me...
    for (int i = beg; i < end; ++i, ++itrack)
    {
      float chi2 = fabs(outChi2[itrack]);//fixme negative chi2 sometimes...
      std::cout << "chi2=" << chi2 << " minChi2[itrack]=" << minChi2[itrack] << std::endl;      
      if (chi2<minChi2[itrack]) 
      {
	minChi2[itrack]=chi2;
	bestHit[itrack]=ih;
      }
    }

  }//end loop over hits

  //copy in MkFitter the hit with lowest chi2
  int itrack = 0;
  //fixme: please vectorize me...
  for (int i = beg; i < end; ++i, ++itrack)
  {
    //fixme decide what to do in case no hit found
    Hit &hit = bestHit[itrack]>=0 ? lay_hits[ bestHit[itrack] ] : lay_hits[0];
    float& chi2 = bestHit[itrack]>=0 ? minChi2[itrack] : minChi2[0];

    std::cout << "ADD BEST HIT FOR TRACK #" << i << std::endl;
    std::cout << "prop x=" << Par[iP].ConstAt(itrack, 0, 0) << " y=" << Par[iP].ConstAt(itrack, 1, 0) << std::endl;      
    std::cout << "copy in hit #" << bestHit[itrack] << " x=" << hit.position()[0] << " y=" << hit.position()[1] << std::endl;    

    msErr[Nhits].CopyIn(itrack, hit.error().Array());
    msPar[Nhits].CopyIn(itrack, hit.parameters().Array());
    Chi2(itrack, 0, 0) += chi2;
  }

  //now update the track parameters with this hit (note that some calculations are already done when computing chi2... not sure it's worth caching them?)
  std::cout << "update parameters" << std::endl;
  updateParametersMPlex(Err[iP], Par[iP], msErr[Nhits], msPar[Nhits],
			Err[iC], Par[iC]);

}
