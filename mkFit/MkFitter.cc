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

void MkFitter::InputTracksAndHitIdx(std::vector<std::vector<Track> >& tracks, std::vector<std::pair<int,int> >& idxs, int beg, int end)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {

    Track &trk = tracks[idxs[i].first][idxs[i].second];

    SeedIdx(itrack, 0, 0) = idxs[i].first;

    Err[iC].CopyIn(itrack, trk.errors().Array());
    Par[iC].CopyIn(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();

    for (int hi = 0; hi < Nhits; ++hi)
    {

      HitsIdx[hi](itrack, 0, 0) = trk.getHitIdx(hi);//dummy value for now

    }
  }
}

int MkFitter::countInvalidHits(int itrack)
{
  int result = 0;
  for (int hi = 0; hi < Nhits; ++hi)
    {
      if (HitsIdx[hi](itrack, 0, 0)<0) result++;
    }
  return result;
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
      tracks[i].addHitIdx(0,0.);//fixme

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

    Hit &hit = lay_hits[ih];
#ifdef DEBUG
    std::cout << "consider hit #" << ih << " of " << lay_hits.size() << std::endl;
    std::cout << "hit x=" << hit.position()[0] << " y=" << hit.position()[1] << std::endl;      
#endif

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
#ifdef DEBUG
      std::cout << "chi2=" << chi2 << " minChi2[itrack]=" << minChi2[itrack] << std::endl;      
#endif
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

#ifdef DEBUG
    std::cout << "ADD BEST HIT FOR TRACK #" << i << std::endl;
    std::cout << "prop x=" << Par[iP].ConstAt(itrack, 0, 0) << " y=" << Par[iP].ConstAt(itrack, 1, 0) << std::endl;      
    std::cout << "copy in hit #" << bestHit[itrack] << " x=" << hit.position()[0] << " y=" << hit.position()[1] << std::endl;    
#endif

    msErr[Nhits].CopyIn(itrack, hit.error().Array());
    msPar[Nhits].CopyIn(itrack, hit.parameters().Array());
    Chi2(itrack, 0, 0) += chi2;
  }

  //now update the track parameters with this hit (note that some calculations are already done when computing chi2... not sure it's worth caching them?)
#ifdef DEBUG
  std::cout << "update parameters" << std::endl;
#endif
  updateParametersMPlex(Err[iP], Par[iP], msErr[Nhits], msPar[Nhits],
			Err[iC], Par[iC]);

}

void MkFitter::FindCandidates(std::vector<Hit>& lay_hits, int firstHit, int lastHit, int beg, int end, std::vector<std::vector<Track> >& tmp_candidates, int offset)
{

  float maxChi2Cut = 30./2.;

  //outer loop over hits, so that tracks can be vectorized
  int ih = firstHit;
  for ( ; ih<lastHit; ++ih)
  {

    Hit &hit = lay_hits[ih];
#ifdef DEBUG
    std::cout << "consider hit #" << ih << " of " << lay_hits.size() << std::endl;
    std::cout << "hit x=" << hit.position()[0] << " y=" << hit.position()[1] << std::endl;      
#endif

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

    //now update the track parameters with this hit (note that some calculations are already done when computing chi2, to be optimized)
    //this is not needed for candidates the hit is not added to, but it's vectorized so doing it serially below should take the same time
    //still it's a waste of time in case the hit is not added to any of the candidates, so check beforehand that at least one cand needs update
    bool oneCandPassCut = false;
    itrack = 0;
    for (int i = beg; i < end; ++i, ++itrack)
      {
	float chi2 = fabs(outChi2[itrack]);//fixme negative chi2 sometimes...
#ifdef DEBUG
	std::cout << "chi2=" << chi2 << std::endl;
#endif
	if (chi2<maxChi2Cut)
	  {
	    oneCandPassCut = true;
	    break;
	  }
      }

    if (oneCandPassCut) { 

      updateParametersMPlex(Err[iP], Par[iP], msErr_oneHit, msPar_oneHit,
			    Err[iC], Par[iC]);
#ifdef DEBUG
      std::cout << "update parameters" << std::endl;
      std::cout << "propagated track parameters x=" << Par[iP].ConstAt(itrack, 0, 0) << " y=" << Par[iP].ConstAt(itrack, 1, 0) << std::endl;
      std::cout << "               hit position x=" << msPar[iP].ConstAt(itrack, 0, 0) << " y=" << msPar[iP].ConstAt(itrack, 1, 0) << std::endl;
      std::cout << "   updated track parameters x=" << Par[iC].ConstAt(itrack, 0, 0) << " y=" << Par[iC].ConstAt(itrack, 1, 0) << std::endl;
#endif

      //create candidate with hit in case chi2<maxChi2Cut
      itrack = 0;
      //fixme: please vectorize me... (not sure it's possible in this case)
      for (int i = beg; i < end; ++i, ++itrack)
	{
	  float chi2 = fabs(outChi2[itrack]);//fixme negative chi2 sometimes...
#ifdef DEBUG
	  std::cout << "chi2=" << chi2 << std::endl;      
#endif
	  if (chi2<maxChi2Cut)
	    {
#ifdef DEBUG
	      std::cout << "chi2 cut passed, creating new candidate" << std::endl;
#endif
	      //create a new candidate and fill the reccands_tmp vector
	      Track newcand;
	      newcand.resetHits();//probably not needed
	      newcand.charge() = Chg(itrack, 0, 0);
	      newcand.setChi2(Chi2(itrack, 0, 0));
	      for (int hi = 0; hi < Nhits; ++hi)
		{
		  newcand.addHitIdx(HitsIdx[hi](itrack, 0, 0),0.);//this should be ok since we already set the chi2 above
		}
#ifdef DEBUG
	      std::cout << "output new hit with x=" << hit.position()[0] << std::endl;
#endif

	      newcand.addHitIdx(ih,chi2);
	      //set the track state to the updated parameters
	      Err[iC].CopyOut(itrack, newcand.errors().Array());
	      Par[iC].CopyOut(itrack, newcand.parameters().Array());
	      
#ifdef DEBUG
	      std::cout << "updated track parameters x=" << newcand.parameters()[0] << " y=" << newcand.parameters()[1] << std::endl;
#endif
	      
	      tmp_candidates[SeedIdx(itrack, 0, 0)-offset].push_back(newcand);
	    }
	}
    }//end if (oneCandPassCut)

  }//end loop over hits

  //now add invalid hit
  int itrack = 0;
  //fixme: please vectorize me...
  for (int i = beg; i < end; ++i, ++itrack)
    {
      if (countInvalidHits(itrack)>0) continue;//check this is ok for vectorization //fixme not optimal
      Track newcand;
      newcand.resetHits();//probably not needed
      newcand.charge() = Chg(itrack, 0, 0);
      newcand.setChi2(Chi2(itrack, 0, 0));
      for (int hi = 0; hi < Nhits; ++hi)
	{
	  newcand.addHitIdx(HitsIdx[hi](itrack, 0, 0),0.);//this should be ok since we already set the chi2 above
	}
      newcand.addHitIdx(-1,0.);
      //set the track state to the propagated parameters
      Err[iP].CopyOut(itrack, newcand.errors().Array());
      Par[iP].CopyOut(itrack, newcand.parameters().Array());	      
      tmp_candidates[SeedIdx(itrack, 0, 0)-offset].push_back(newcand);
    }

}


void MkFitter::GetHitRange(std::vector<BinInfo>& segmentMapLay_, int beg, int end, const float& etaDet, int& firstHit, int& lastHit) {

    int itrack = 0;
    //cannot be vectorized, I think
    for (int i = beg; i < end; ++i, ++itrack)
      {
	float eta = getEta(Par[iP].ConstAt(itrack, 3, 0),Par[iP].ConstAt(itrack, 4, 0),Par[iP].ConstAt(itrack, 5, 0));
	//protect against anomalous eta (should go into getEtaPartition maybe?)
	if (fabs(eta)>etaDet) eta = (eta>0 ? etaDet*0.99 : -etaDet*0.99);
	unsigned int etabin = getEtaPartition(eta,etaDet);
	BinInfo binInfo = segmentMapLay_[etabin];
#ifdef DEBUG
	std::cout << "propagated track parameters eta=" << eta << " bin=" << etabin << " begin=" << binInfo.first << " size=" << binInfo.second << std::endl;
#endif
	if (firstHit==-1 || binInfo.first<firstHit) firstHit =  binInfo.first;
	if (lastHit==-1 || (binInfo.first+binInfo.second)>lastHit) lastHit = binInfo.first+binInfo.second;
      }
#ifdef DEBUG
    std::cout << "found range firstHit=" << firstHit << " lastHit=" << lastHit << std::endl;
#endif
}
