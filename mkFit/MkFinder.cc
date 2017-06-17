#include "MkFinder.h"

#include "CandCloner.h"
#include "SteeringParams.h"
#include "HitStructures.h"

#include "KalmanUtilsMPlex.h"

//#define DEBUG
#include "Debug.h"

//==============================================================================
// Input / Output TracksAndHitIdx
//==============================================================================

void MkFinder::InputTracksAndHitIdx(const std::vector<Track>& tracks,
                                    int beg, int end, bool inputProp)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  const int iI = inputProp ? iP : iC;

  for (int i = beg, imp = 0; i < end; ++i, ++imp)
  {
    copy_in(tracks[i], imp, iI);
  }
}

void MkFinder::InputTracksAndHitIdx(const std::vector<Track>& tracks,
                                    const std::vector<int>  & idxs,
                                    int beg, int end, bool inputProp, int mp_offset)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  const int iI = inputProp ? iP : iC;

  for (int i = beg, imp = mp_offset; i < end; ++i, ++imp)
  {
    copy_in(tracks[idxs[i]], imp, iI);
  }
}

void MkFinder::InputTracksAndHitIdx(const std::vector<std::vector<Track> >& tracks,
                                    const std::vector<std::pair<int,int> >& idxs,
                                    int beg, int end, bool inputProp)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  const int iI = inputProp ? iP : iC;

  for (int i = beg, imp = 0; i < end; ++i, ++imp)
  {
    const Track &trk = tracks[idxs[i].first][idxs[i].second];

    copy_in(trk, imp, iI);

    SeedIdx(imp, 0, 0) = idxs[i].first;
    CandIdx(imp, 0, 0) = idxs[i].second;
  }
}

void MkFinder::InputTracksAndHitIdx(const std::vector<std::vector<Track> >& tracks,
                                    const std::vector<std::pair<int,MkFinder::IdxChi2List> >& idxs,
                                    int beg, int end, bool inputProp)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  const int iI = inputProp ? iP : iC;

  for (int i = beg, imp = 0; i < end; ++i, ++imp)
  {
    const Track &trk = tracks[idxs[i].first][idxs[i].second.trkIdx];

    copy_in(trk, imp, iI);

    SeedIdx(imp, 0, 0) = idxs[i].first;
    CandIdx(imp, 0, 0) = idxs[i].second.trkIdx;
  }
}

void MkFinder::OutputTracksAndHitIdx(std::vector<Track>& tracks,
                                     int beg, int end, bool outputProp) const
{
  // Copies requested track parameters into Track objects.
  // The tracks vector should be resized to allow direct copying.

  const int iO = outputProp ? iP : iC;

  for (int i = beg, imp = 0; i < end; ++i, ++imp)
  {
    copy_out(tracks[i], imp, iO);
  }
}

void MkFinder::OutputTracksAndHitIdx(std::vector<Track>& tracks,
                                     const std::vector<int>& idxs,
                                     int beg, int end, bool outputProp) const
{
  // Copies requested track parameters into Track objects.
  // The tracks vector should be resized to allow direct copying.

  const int iO = outputProp ? iP : iC;

  for (int i = beg, imp = 0; i < end; ++i, ++imp)
  {
    copy_out(tracks[idxs[i]], imp, iO);
  }
}


//==============================================================================
// SelectHitIndices
//==============================================================================

void MkFinder::SelectHitIndices(const LayerOfHits &layer_of_hits,
                                const int N_proc, bool dump)
{
  // debug = 1;
  // dump = true;

  const LayerOfHits &L = layer_of_hits;
  const int   iI = iP;
  const float nSigmaPhi = 3;
  const float nSigmaZ   = 3;
  const float nSigmaR   = 3;

  if (dump)
    printf("LayerOfHits::SelectHitIndices %s layer=%d N_proc=%d\n",
           L.is_barrel() ? "barrel" : "endcap", L.layer_id(), N_proc);

  // Vectorizing this makes it run slower!
  //#pragma ivdep
  //#pragma simd
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    XHitSize[itrack] = 0;

    const float x = Par[iI].ConstAt(itrack, 0, 0);
    const float y = Par[iI].ConstAt(itrack, 1, 0);

    const float r2     = x*x + y*y;
    const float dphidx = -y/r2, dphidy = x/r2;
    const float dphi2  = dphidx * dphidx * Err[iI].ConstAt(itrack, 0, 0) +
      dphidy * dphidy * Err[iI].ConstAt(itrack, 1, 1) +
      2 * dphidx * dphidy * Err[iI].ConstAt(itrack, 0, 1);
#ifdef HARD_CHECK
    assert(dphi2 >= 0);
#endif

    float q, dq, phi, dphi;

    phi  = getPhi(x, y);
    dphi = nSigmaPhi * std::sqrt(std::abs(dphi2));
    dphi = std::max(std::abs(dphi), L.min_dphi());

    if (L.is_barrel())
    {
      float z  = Par[iI].ConstAt(itrack, 2, 0);
      float dz = nSigmaZ * std::sqrt(Err[iI].ConstAt(itrack, 2, 2));

      dz = std::max(std::abs(dz),   L.min_dq());

      if (Config::useCMSGeom)
      {
        //now correct for bending and for layer thickness unsing linear approximation
        //fixme! using constant value, to be taken from layer properties
        //XXXXMT4GC should we also increase dz?
        const float deltaR = Config::cmsDeltaRad;
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

      q =  z;
      dq = dz;
    }
    else // endcap
    {
      float  r = std::sqrt(r2);
      float dr = nSigmaR*(x*x*Err[iI].ConstAt(itrack, 0, 0) + y*y*Err[iI].ConstAt(itrack, 1, 1) + 2*x*y*Err[iI].ConstAt(itrack, 0, 1))/r2;

      dr = std::max(std::abs(dr), L.min_dq());

      if (Config::useCMSGeom)
      {
#ifdef CCSCOORD
        //now correct for bending and for layer thickness unsing linear approximation
        //fixme! using constant value, to be taken from layer properties
        //XXXXMT4GC should we also increase dr?
        const float deltaZ = 5;
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

      q =  r;
      dq = dr;
    }

    dphi = std::min(std::abs(dphi), L.max_dphi());
    dq   = std::min(std::abs(dq),   L.max_dq());

    const int qb1 = L.GetQBinChecked(q - dq);
    const int qb2 = L.GetQBinChecked(q + dq) + 1;
    const int pb1 = L.GetPhiBin(phi - dphi);
    const int pb2 = L.GetPhiBin(phi + dphi) + 1;
    // MT: The extra phi bins give us ~1.5% more good tracks at expense of 10% runtime.
    // const int pb1 = L.GetPhiBin(phi - dphi) - 1;
    // const int pb2 = L.GetPhiBin(phi + dphi) + 2;

    if (dump)
      printf("  %2d: %6.3f %6.3f %6.6f %7.5f %3d %3d %4d %4d\n",
             itrack, q, phi, dq, dphi, qb1, qb2, pb1, pb2);

    // MT: One could iterate in "spiral" order, to pick hits close to the center.
    // http://stackoverflow.com/questions/398299/looping-in-a-spiral
    // This would then work best with relatively small bin sizes.
    // Or, set them up so I can always take 3x3 array around the intersection.

    for (int qi = qb1; qi < qb2; ++qi)
    {
      for (int pi = pb1; pi < pb2; ++pi)
      {
        const int pb = pi & L.m_phi_mask;

        // MT: The following line is the biggest hog (4% total run time).
        // This comes from cache misses, I presume.
        // It might make sense to make first loop to extract bin indices
        // and issue prefetches at the same time.
        // Then enter vectorized loop to actually collect the hits in proper order.

        for (int hi = L.m_phi_bin_infos[qi][pb].first; hi < L.m_phi_bin_infos[qi][pb].second; ++hi)
        {
          // MT: Access into m_hit_zs and m_hit_phis is 1% run-time each.

#ifdef LOH_USE_PHI_Q_ARRAYS
          float ddq   = std::abs(q   - L.m_hit_qs[hi]);
          float ddphi = std::abs(phi - L.m_hit_phis[hi]);
          if (ddphi > Config::PI) ddphi = Config::TwoPI - ddphi;

          if (dump)
            printf("     SHI %3d %4d %4d %5d  %6.3f %6.3f %6.4f %7.5f   %s\n",
                   qi, pi, pb, hi,
                   L.m_hit_qs[hi], L.m_hit_phis[hi], ddq, ddphi,
                   (ddq < dq && ddphi < dphi) ? "PASS" : "FAIL");

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
// AddBestHit - Best Hit Track Finding
//==============================================================================

//#define NO_PREFETCH
//#define NO_GATHER

void MkFinder::AddBestHit(const LayerOfHits    &layer_of_hits, const int N_proc,
                          const SteeringParams &st_par)
{
  // debug = true;

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
        msErr.CopyIn(itrack, hit.errArray());
        msPar.CopyIn(itrack, hit.posArray());
      }
    }

#else //NO_GATHER

#if defined(MIC_INTRINSICS)
    msErr.SlurpIn(varr + off_error, vi);
    msPar.SlurpIn(varr + off_param, vi);
#else
    msErr.SlurpIn(varr + off_error, idx);
    msPar.SlurpIn(varr + off_param, idx);
#endif
#endif //NO_GATHER

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    (*st_par.compute_chi2_foo)(Err[iP], Par[iP], Chg, msErr, msPar,
                               outChi2, N_proc);

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
      const Hit  &hit  = layer_of_hits.m_hits[ bestHit[itrack] ];
      const float chi2 = minChi2[itrack];

      dprint("ADD BEST HIT FOR TRACK #" << itrack << std::endl
        << "prop x=" << Par[iP].ConstAt(itrack, 0, 0) << " y=" << Par[iP].ConstAt(itrack, 1, 0) << std::endl
        << "copy in hit #" << bestHit[itrack] << " x=" << hit.position()[0] << " y=" << hit.position()[1]);

      msErr.CopyIn(itrack, hit.errArray());
      msPar.CopyIn(itrack, hit.posArray());
      Chi2(itrack, 0, 0) += chi2;

      add_hit(itrack, bestHit[itrack], layer_of_hits.layer_id());
    }
    else
    {
      int fake_hit_idx = -1;

      if (layer_of_hits.is_endcap() &&
          layer_of_hits.is_in_xy_hole(Par[iP](itrack,0,0), Par[iP](itrack,1,0)))
      {
        // YYYYYY Config::store_missed_layers
        fake_hit_idx = -3;
      }

      dprint("ADD FAKE HIT FOR TRACK #" << itrack << " withinBounds=" << (fake_hit_idx != -3) << " r=" << std::hypot(Par[iP](itrack,0,0), Par[iP](itrack,1,0)));

      msErr.SetDiagonal3x3(itrack, 666);
      msPar(itrack,0,0) = Par[iP](itrack,0,0);
      msPar(itrack,1,0) = Par[iP](itrack,1,0);
      msPar(itrack,2,0) = Par[iP](itrack,2,0);
      // Don't update chi2

      add_hit(itrack, fake_hit_idx, layer_of_hits.layer_id());
    }
  }

  // Update the track parameters with this hit. (Note that some calculations
  // are already done when computing chi2. Not sure it's worth caching them?)

  dprint("update parameters");
  (*st_par.update_param_foo)(Err[iP], Par[iP], Chg, msErr, msPar,
                             Err[iC], Par[iC], N_proc);

  //std::cout << "Par[iP](0,0,0)=" << Par[iP](0,0,0) << " Par[iC](0,0,0)=" << Par[iC](0,0,0)<< std::endl;
}

//==============================================================================
// FindCandidates - Standard Track Finding
//==============================================================================

void MkFinder::FindCandidates(const LayerOfHits &layer_of_hits,
                              std::vector<std::vector<Track>>& tmp_candidates,
                              const int offset, const int N_proc,
                              const SteeringParams &st_par)
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
    msErr.SlurpIn(varr + off_error, vi);
    msPar.SlurpIn(varr + off_param, vi);
#else
    msErr.SlurpIn(varr + off_error, idx);
    msPar.SlurpIn(varr + off_param, idx);
#endif

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    (*st_par.compute_chi2_foo)(Err[iP], Par[iP], Chg, msErr, msPar,
                               outChi2, N_proc);

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
      (*st_par.update_param_foo)(Err[iP], Par[iP], Chg, msErr, msPar,
                                 Err[iC], Par[iC], N_proc);

      dprint("update parameters" << std::endl
	     << "propagated track parameters x=" << Par[iP].ConstAt(0, 0, 0) << " y=" << Par[iP].ConstAt(0, 1, 0) << std::endl
	     << "               hit position x=" << msPar.ConstAt(0, 0, 0)   << " y=" << msPar.ConstAt(0, 1, 0) << std::endl
	     << "   updated track parameters x=" << Par[iC].ConstAt(0, 0, 0) << " y=" << Par[iC].ConstAt(0, 1, 0));

      //create candidate with hit in case chi2 < Config::chi2Cut
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
            copy_out(newcand, itrack, iC);
	    newcand.addHitIdx(XHitArr.At(itrack, hit_cnt, 0), layer_of_hits.layer_id(), chi2);

	    dprint("updated track parameters x=" << newcand.parameters()[0] << " y=" << newcand.parameters()[1] << " z=" << newcand.parameters()[2] << " pt=" << 1./newcand.parameters()[3]);

	    tmp_candidates[SeedIdx(itrack, 0, 0) - offset].emplace_back(newcand);
	  }
	}
      }
    }//end if (oneCandPassCut)

  }//end loop over hits

  //now add invalid hit
  //fixme: please vectorize me...
  for (int itrack = 0; itrack < N_proc; ++itrack)
  {
    int fake_hit_idx = num_invalid_hits(itrack) < Config::maxHolesPerCand ? -1 : -2;

    if (layer_of_hits.is_endcap() &&
        layer_of_hits.is_in_xy_hole(Par[iP](itrack,0,0), Par[iP](itrack,1,0)))
    {
      // YYYYYY Config::store_missed_layers
      fake_hit_idx = -3;
    }

    dprint("ADD FAKE HIT FOR TRACK #" << itrack << " withinBounds=" << (fake_hit_idx != -3) << " r=" << std::hypot(Par[iP](itrack,0,0), Par[iP](itrack,1,0)));

    Track newcand;
    copy_out(newcand, itrack, iP);
    newcand.addHitIdx(fake_hit_idx, layer_of_hits.layer_id(), 0.);
    tmp_candidates[SeedIdx(itrack, 0, 0) - offset].emplace_back(newcand);
  }
}


//==============================================================================
// FindCandidatesCloneEngine - Clone Engine Track Finding
//==============================================================================

void MkFinder::FindCandidatesCloneEngine(const LayerOfHits &layer_of_hits, CandCloner& cloner,
                                         const int offset, const int N_proc,
                                         const SteeringParams &st_par)
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
    msErr.SlurpIn(varr + off_error, vi);
    msPar.SlurpIn(varr + off_param, vi);
#else
    msErr.SlurpIn(varr + off_error, idx);
    msPar.SlurpIn(varr + off_param, idx);
#endif

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    (*st_par.compute_chi2_foo)(Err[iP], Par[iP], Chg, msErr, msPar, outChi2, N_proc);

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
          tmpList.nhits  = NFoundHits(itrack,0,0) + 1;
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

    int fake_hit_idx = num_invalid_hits(itrack) < Config::maxHolesPerCand ? -1 : -2;

    if (layer_of_hits.is_endcap() &&
          layer_of_hits.is_in_xy_hole(Par[iP](itrack,0,0), Par[iP](itrack,1,0)))
    {
      // YYYYYY Config::store_missed_layers
      fake_hit_idx = -3;
    }

    IdxChi2List tmpList;
    tmpList.trkIdx = CandIdx(itrack, 0, 0);
    tmpList.hitIdx = fake_hit_idx;
    tmpList.nhits  = NFoundHits(itrack,0,0);
    tmpList.chi2   = Chi2(itrack, 0, 0);
    cloner.add_cand(SeedIdx(itrack, 0, 0) - offset, tmpList);
#ifdef DEBUG
    std::cout << "adding invalid hit" << std::endl;
#endif
  }
}


//==============================================================================
// UpdateWithLastHit
//==============================================================================

void MkFinder::UpdateWithLastHit(const LayerOfHits &layer_of_hits, int N_proc,
                                 const SteeringParams &st_par)
{
  // layer_of_hits is for previous layer!

  for (int i = 0; i < N_proc; ++i)
  {
    int hit_idx = HoTArrs[i][ NHits[i] - 1].index;

    if (hit_idx < 0) continue;

    Hit &hit = layer_of_hits.m_hits[hit_idx];

    msErr.CopyIn(i, hit.errArray());
    msPar.CopyIn(i, hit.posArray());
  }

  (*st_par.update_param_foo)(Err[iP], Par[iP], Chg, msErr, msPar,
                             Err[iC], Par[iC], N_proc);

  //now that we have moved propagation at the end of the sequence we lost the handle of
  //using the propagated parameters instead of the updated for the missing hit case.
  //so we need to replace by hand the updated with the propagated
  //there may be a better way to restore this...

  for (int i = 0; i < N_proc; ++i)
  {
    if (HoTArrs[i][ NHits[i] - 1].index < 0)
    {
      float tmp[21];
      Err[iP].CopyOut(i, tmp);
      Err[iC].CopyIn (i, tmp);
      Par[iP].CopyOut(i, tmp);
      Par[iC].CopyIn (i, tmp);
    }
  }
}


//==============================================================================
// CopyOutParErr
//==============================================================================

void MkFinder::CopyOutParErr(std::vector<std::vector<Track> >& seed_cand_vec,
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
              << " posEta=" << cand.posEta());
  }
}
