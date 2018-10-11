#include "MkFinderFV.h"

#include "SteeringParams.h"
#include "HitStructures.h"

#include "KalmanUtilsMPlex.h"
#include "MatriplexPackers.h"

#include <set>

//#define DEBUG
#include "Debug.h"

namespace mkfit {

struct candComp {
  bool operator()(const MkFinderFv::IdxChi2List& c1, const MkFinderFv::IdxChi2List& c2) const
  {
    if (c1.nhits == c2.nhits) return c1.chi2 < c2.chi2;
    return c1.nhits > c2.nhits;
  }
};

// fill all candidates with the input seed
template<int nseeds, int ncands>
void MkFinderFV<nseeds, ncands>::InputTrack(const Track& track, int iseed, int offset, bool inputProp)
{
  const int iI = inputProp ? iP : iC;
  for (int i = 0; i < ncands; ++i)
  {
    int imp = i + offset*ncands;
    copy_in(track, imp, iI);
    SeedIdx(imp, 0, 0) = iseed;
    CandIdx(imp, 0, 0) = i == 0 ? 0 : -1; // dummies are -1
  }
}

template<int nseeds, int ncands>
void MkFinderFV<nseeds, ncands>::OutputTrack(std::vector<Track>& tracks,
                           int itrack, int imp, bool outputProp) const
{
  // Copies requested track parameters into Track object.
  // The tracks vector should be resized to allow direct copying.

  const int iO = outputProp ? iP : iC;
  copy_out(tracks[itrack], imp, iO);
}

//==============================================================================
// SelectHitIndices
//==============================================================================

template<int nseeds, int ncands>
void MkFinderFV<nseeds, ncands>::SelectHitIndices(const LayerOfHits &layer_of_hits)
{
  //debug = true;

  const LayerOfHits &L = layer_of_hits;
  const int   iI = iP;
  const float nSigmaPhi = 3;
  const float nSigmaZ   = 3;
  const float nSigmaR   = 3;

  dprintf("LayerOfHits::SelectHitIndices %s layer=%d\n",
           L.is_barrel() ? "barrel" : "endcap", L.layer_id());

  float dqv[NNFV], dphiv[NNFV], qv[NNFV], phiv[NNFV];
  int goodv[NNFV];

  int best[nseeds];
  for (auto i = 0; i < nseeds; ++i) { best[i] = std::max(BestCandidate(i), 0); }

  const auto assignbins = [&](int itrack, float q, float dq, float phi, float dphi){
    dphi = std::min(std::abs(dphi), L.max_dphi());
    dq   = clamp(dq, L.min_dq(), L.max_dq());

    dprintf("     SHI %5d  %6.3f %6.3f %6.4f %7.5f\n", itrack, q, phi, dq, dphi);

    qv[itrack] = q;
    phiv[itrack] = phi;
    dphiv[itrack] = dphi;
    dqv[itrack]   = dq;
    goodv[itrack] = XWsrResult[itrack].m_wsr != WSR_Outside && CandIdx[itrack] >= 0;
    // MT: The extra phi bins give us ~1.5% more good tracks at expense of 10% runtime.
    //     That was for 10k cylindrical cow, I think.
    // const int pb1 = L.GetPhiBin(phi - dphi) - 1;
    // const int pb2 = L.GetPhiBin(phi + dphi) + 2;
  };

  const auto calcdphi2 = [&](int itrack, float dphidx, float dphidy) {
    return dphidx * dphidx * Err[iI].ConstAt(itrack, 0, 0) +
           dphidy * dphidy * Err[iI].ConstAt(itrack, 1, 1) +
       2 * dphidx * dphidy * Err[iI].ConstAt(itrack, 0, 1);
  };

  const auto calcdphi = [&](float dphi2) {
    return std::max(nSigmaPhi * std::sqrt(std::abs(dphi2)), L.min_dphi());
  };

  // Pull out the part of the loop that vectorizes
  if (L.is_barrel())
  {
#pragma omp simd
    for (int itrack = 0; itrack < NNFV; ++itrack)
    {
      XHitSize[itrack] = 0;

      const float x = Par[iI].ConstAt(itrack, 0, 0);
      const float y = Par[iI].ConstAt(itrack, 1, 0);

      const float r2     = x*x + y*y;
      const float dphidx = -y/r2, dphidy = x/r2;
      const float dphi2  = calcdphi2(itrack, dphidx, dphidy);
  #ifdef HARD_CHECK
      assert(dphi2 >= 0);
  #endif

      const float phi  = getPhi(x, y);
      float dphi = calcdphi(dphi2);

      const float z  = Par[iI].ConstAt(itrack, 2, 0);
      const float dz = std::abs(nSigmaZ * std::sqrt(Err[iI].ConstAt(itrack, 2, 2)));

      if (Config::useCMSGeom)
      {
        //now correct for bending and for layer thickness unsing linear approximation
        //fixme! using constant value, to be taken from layer properties
        //XXXXMT4GC should we also increase dz?
        const float deltaR = Config::cmsDeltaRad;
        const float r  = std::sqrt(r2);
        //here alpha is the difference between posPhi and momPhi
        const float alpha = phi - Par[iP].ConstAt(itrack, 4, 0);
        float cosA, sinA;
        if (Config::useTrigApprox) {
          sincos4(alpha, sinA, cosA);
        } else {
          cosA = std::cos(alpha);
          sinA = std::sin(alpha);
        }
        //take abs so that we always inflate the window
        const float dist = std::abs(deltaR*sinA/cosA);
        dphi += dist / r;
      }

      XWsrResult[itrack] = L.is_within_z_sensitive_region(z, dz);
      assignbins(itrack, z, dz, phi, dphi);
    }
  } else {
#pragma omp simd
    for (int itrack = 0; itrack < NNFV; ++itrack)
    {
      XHitSize[itrack] = 0;

      const float x = Par[iI].ConstAt(itrack, 0, 0);
      const float y = Par[iI].ConstAt(itrack, 1, 0);

      const float r2     = x*x + y*y;
      const float dphidx = -y/r2, dphidy = x/r2;
      const float dphi2  = calcdphi2(itrack, dphidx, dphidy);
  #ifdef HARD_CHECK
      assert(dphi2 >= 0);
  #endif

      const float phi  = getPhi(x, y);
      float dphi = calcdphi(dphi2);

      const float  r = std::sqrt(r2);
      const float dr = std::abs(nSigmaR*(x*x*Err[iI].ConstAt(itrack, 0, 0) + y*y*Err[iI].ConstAt(itrack, 1, 1) + 2*x*y*Err[iI].ConstAt(itrack, 0, 1)) / r2);

      if (Config::useCMSGeom)
      {
        //now correct for bending and for layer thickness unsing linear approximation
        //fixme! using constant value, to be taken from layer properties
        //XXXXMT4GC should we also increase dr?
        const float deltaZ = 5;
        const float cosT = std::cos(Par[iI].ConstAt(itrack, 5, 0));
        const float sinT = std::sin(Par[iI].ConstAt(itrack, 5, 0));
        //here alpha is the helix angular path corresponding to deltaZ
        const float k = Chg.ConstAt(itrack, 0, 0) * 100.f / (-Config::sol*Config::Bfield);
        const float alpha  = deltaZ*sinT*Par[iI].ConstAt(itrack, 3, 0)/(cosT*k);
        dphi += std::abs(alpha);
      }

      XWsrResult[itrack] = L.is_within_r_sensitive_region(r, dr);
      assignbins(itrack, r, dr, phi, dphi);
    }
  }

  int qb1[nseeds], qb2[nseeds], pb1[nseeds], pb2[nseeds];
  for (auto iseed = 0; iseed < nseeds; ++iseed) {
    const float qbest   =   qv[best[iseed]];
    const float phibest = phiv[best[iseed]];
    float dqvar = 0.0f, dphivar = 0.0f;
    int count = 0;
    const int base = iseed*ncands;
    #pragma omp simd
    for (int i = 0; i < ncands; ++i) {
      const int itrack = base + i;
      const int weight = goodv[itrack];
      count   += weight;
      dqvar   += weight*std::abs(  qv[itrack] - qbest);
      dphivar += weight*std::abs(phiv[itrack] - phibest);
    }
    float c = 1.0f;
    if (count > 1) {
      c = 1.0f/(count-1);
    }
    const float dq   =   c*dqvar + dqv[best[iseed]];
    const float dphi = c*dphivar + dphiv[best[iseed]];

    qb1[iseed] = L.GetQBinChecked(qbest - dq);
    qb2[iseed] = L.GetQBinChecked(qbest + dq) + 1;
    pb1[iseed] = L.GetPhiBin(phibest - dphi);
    pb2[iseed] = L.GetPhiBin(phibest + dphi) + 1;
    dprint("Seed " << iseed << " Best " << BestCandidate(iseed) << " q " << qbest << " phi " << phibest 
      << " c " << c << " dqvar " << dqvar << " dphivar " << dphivar << " bins " 
      << qb1[iseed] << "-" << qb2[iseed] << ", "
      << pb1[iseed] << "-" << pb2[iseed]);
  }

  _mm_prefetch((const char*) &L.m_phi_bin_infos[qb1[0]][pb1[0] & L.m_phi_mask], _MM_HINT_T0);

  for (auto iseed = 0; iseed < nseeds; ++iseed) {
    const int base = iseed*ncands;
    for (int qi = qb1[iseed]; qi < qb2[iseed]; ++qi) {
      if (qi+1 < qb2[iseed]) {
        _mm_prefetch((const char*) &L.m_phi_bin_infos[qi+1][pb1[iseed] & L.m_phi_mask], _MM_HINT_T0);
      }
      for (int pi = pb1[iseed]; pi < pb2[iseed]; ++pi) {
        int pb = pi & L.m_phi_mask;
        for (int hi = L.m_phi_bin_infos[qi][pb].first; hi < L.m_phi_bin_infos[qi][pb].second; ++hi) {
          unsigned int pass[ncands] = {};
          if (Config::usePhiQArrays) {
            #pragma omp simd
            for (int i = 0; i < ncands; ++i) {
              const int itrack = base + i;
              const float q    = qv[itrack];
              const float dq   = dqv[itrack];
              
              const float phi  = phiv[itrack];
              const float dphi = dphiv[itrack];
              
              const float ddq   =       std::abs(q   - L.m_hit_qs[hi]);
              const float ddphi = cdist(std::abs(phi - L.m_hit_phis[hi]));
              
              dprintf("     SHI %5d  %6.3f %6.3f %6.4f %7.5f   %s\n",
                      hi,
                      L.m_hit_qs[hi], L.m_hit_phis[hi], q, phi,
                      (ddq < dq && ddphi < dphi) ? "PASS" : "FAIL");
              
                      pass[i] = goodv[itrack] && ddq < dq && ddphi < dphi;
            }
          } else {
            #pragma omp simd
            for (int i = 0; i < ncands; ++i) {
              const int itrack = base + i;
              pass[i] = goodv[itrack];
            }
          }
          for (int i = 0; i < ncands; ++i) {
            const int itrack = base + i;
            if (pass[i] && XHitSize[itrack] < MPlexHitIdxMax) {
              XHitArr.At(itrack, XHitSize[itrack]++, 0) = hi;
            }
          }
        }
      }
    }
  }
}


//==============================================================================
// FindCandidates - FV Track Finding
//==============================================================================

template<int nseeds, int ncands>
void MkFinderFV<nseeds, ncands>::FindCandidates(const LayerOfHits &layer_of_hits,
                                                const FindingFoos &fnd_foos)
{
  MatriplexHitPacker mhp;

  const char *varr      = (char*) layer_of_hits.m_hits;

  // prefetch the first set of hits to L1 and the second one to L2.
  #pragma omp simd
  for (int it = 0; it < NNFV; ++it)
  {
    if (XHitSize[it] > 0)
    {
      _mm_prefetch(varr + XHitArr.At(it, 0, 0) * sizeof(Hit), _MM_HINT_T0);
      if (XHitSize[it] > 1)
      {
        _mm_prefetch(varr + XHitArr.At(it, 1, 0) * sizeof(Hit), _MM_HINT_T1);
      }
    }
  }
  // XXXX MT FIXME: use masks to filter out SlurpIns

  // Has basically no effect, it seems.
  //#pragma noprefetch
  const int maxhit = XHitMax();
  for (int hit_cnt = 0; hit_cnt < maxhit; ++hit_cnt)
  {
    mhp.Reset();

#pragma omp simd
    for (int itrack = 0; itrack < NNFV; ++itrack)
    {
      if (hit_cnt < XHitSize[itrack])
      {
        mhp.AddInputAt(itrack, layer_of_hits.m_hits[ XHitArr.At(itrack, hit_cnt, 0) ]);
      }
    }

    // Prefetch to L2 the hits we'll (probably) process after two loops iterations.
    // Ideally this would be initiated before coming here, for whole bunch_of_hits.m_hits vector.
    for (int itrack = 0; itrack < NNFV; ++itrack)
    {
      if (hit_cnt + 2 < XHitSize[itrack])
      {
        _mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+2, 0)*sizeof(Hit), _MM_HINT_T1);
      }
    }

    mhp.Pack(msErr, msPar);

    //now compute the chi2 of track state vs hit
    (*fnd_foos.m_compute_chi2_foo)(Err[iP], Par[iP], Chg, msErr, msPar, XHitChi2[hit_cnt], NNFV, Config::finding_intra_layer_pflags);

    // Prefetch to L1 the hits we'll (probably) process in the next loop iteration.
    for (int itrack = 0; itrack < NNFV; ++itrack)
    {
      if (hit_cnt + 1 < XHitSize[itrack])
      {
        _mm_prefetch(varr + XHitArr.At(itrack, hit_cnt+1, 0)*sizeof(Hit), _MM_HINT_T0);
      }
    }
  }//end loop over hits
}


//==============================================================================
// UpdateWithLastHit
//==============================================================================

template<int nseeds, int ncands>
void MkFinderFV<nseeds, ncands>::UpdateWithLastHit(const LayerOfHits &layer_of_hits, const FindingFoos &fnd_foos)
{
  bool has_hit[NNFV];

  for (int i = 0; i < NNFV; ++i)
  {
    const HitOnTrack &hot = HoTArrs[i][ NHits[i] - 1];
    const int hit_idx = hot.index;

    has_hit[i] = hit_idx >= 0 && hot.layer == layer_of_hits.layer_id();

    if (has_hit[i]) {
      const Hit& hit = layer_of_hits.m_hits[hit_idx];

      msErr.CopyIn(i, hit.errArray());
      msPar.CopyIn(i, hit.posArray());
    }
  }

  (*fnd_foos.m_update_param_foo)(Err[iP], Par[iP], Chg, msErr, msPar,
                                 Err[iC], Par[iC], NNFV, Config::finding_intra_layer_pflags);

  //now that we have moved propagation at the end of the sequence we lost the handle of
  //using the propagated parameters instead of the updated for the missing hit case.
  //so we need to replace by hand the updated with the propagated
  //there may be a better way to restore this...
  for (int i = 0; i < NNFV; ++i)
  {
    if (!has_hit[i]) {
      Err[iC].CopyIn(i, Err[iP], i);
      Par[iC].CopyIn(i, Par[iP], i);
    }
  }
}


//==============================================================================
// SelectBestCandidates
//==============================================================================

template<int nseeds, int ncands>
void MkFinderFV<nseeds, ncands>::SelectBestCandidates(const LayerOfHits &layer_of_hits)
{
  for (auto iseed = 0; iseed < nseeds; ++iseed) {
    std::array<IdxChi2List, ncands*(MPlexHitIdxMax+1)> idxv;
    int idxvsize = 0;

    auto itrack = iseed*ncands;
    for (int it = 0; it < ncands; ++it, ++itrack) {
      if (CandIdx[itrack] >= 0) {
        for (int hit_cnt = 0; hit_cnt < XHitSize[itrack]; ++hit_cnt) {
          const float chi2 = std::abs(XHitChi2[hit_cnt][itrack]);
          if (chi2 < Config::chi2Cut) {
            dprint("Candidate " << itrack << " adding hit " << hit_cnt << " idx " << XHitArr(itrack,hit_cnt,0));
            idxv[idxvsize++] = IdxChi2List(itrack, XHitArr(itrack,hit_cnt,0), NFoundHits[itrack]+1, Chi2[itrack] + chi2);
          }
        }
        int fake_hit_idx = num_invalid_hits(itrack) < Config::maxHolesPerCand ? -1 : -2;

        if (XWsrResult[itrack].m_wsr == WSR_Outside) {
          fake_hit_idx = -4;
        } else {
          if (XWsrResult[itrack].m_wsr == WSR_Edge) {
            fake_hit_idx = -3; // YYYYYY Config::store_missed_layers
          }
        }
        // add a fake hit
        idxv[idxvsize++] = IdxChi2List(itrack, fake_hit_idx, NFoundHits[itrack], Chi2[itrack]);
      }
    }

    if (idxvsize > Config::maxCandsPerSeed) {
      auto b = idxv.begin();
      std::nth_element(b, b+Config::maxCandsPerSeed-1, b+idxvsize, candComp());
    }

    const int nc = std::min(idxvsize, Config::maxCandsPerSeed);
    dprint("nc = " << nc);

    int filled[ncands];
    int tocopy[ncands];
    int tmpidx[ncands];

    for (int i = 0; i < ncands; ++i) {
      filled[i]  = -1;
      tmpidx[i] = -1;
    }

    // identify tracks that can be reused and candidates that need a track to be copied
    int needscopy = 0;
    for (int ic = 0; ic < nc; ++ic) {
      int it = idxv[ic].trkIdx - iseed*ncands;
      if (filled[it] < 0) {
        filled[it] = ic; // slot it gets candidate ic, no copy
        tmpidx[it] = it; // it or ic?
      } else {
        tocopy[needscopy++] = ic; // candidate ic needs a copy
      }
    }

    // map out the copy destinations
    int copyto[ncands];
    int icopy = 0;
    for (int it = 0; it < ncands; ++it) {
      if (filled[it] < 0) {
        if (icopy < needscopy) { // empty slot
          copyto[icopy++] = it;
        } else {
          filled[it] = 0; // dummy
        }
      }
    }

    // copy the ones that need copying
    for (int i = 0; i < icopy; ++i) {
      int it = copyto[i] + iseed*ncands; // fill this slot

      int ic = tocopy[i];
      int ot = idxv[ic].trkIdx; // filled from this slot
      dprint("copying slot " << ot << " to slot " << it);

      filled[copyto[i]] = ic;

      Err[iP].Copy(it, ot);
      Par[iP].Copy(it, ot);

      NHits     (it, 0, 0) = NHits(ot, 0, 0);
      NFoundHits(it, 0, 0) = NFoundHits(ot, 0, 0);

      SeedIdx[it] = SeedIdx[ot];
      tmpidx[copyto[i]] = it; // it or ic?

      memcpy(HoTArrs[it], HoTArrs[ot], sizeof(HoTArrs[it]));
    }

    // update all
    for (int ic = 0; ic < ncands; ++ic) {
      int it = ic + iseed*ncands;
      CandIdx[it] = tmpidx[ic];
      if (CandIdx[it] >= 0) {
        const auto& idx = idxv[filled[ic]];
        add_hit(it, idx.hitIdx, layer_of_hits.layer_id());
        Chi2(it, 0, 0) = idx.chi2;
#ifdef DEBUG
        if (debug) {
          dmutex_guard;
          std::cout << "slot " << it << " adding hit idx " << std::setw(4) << idx.hitIdx << ":";
          for (int i = 0; i < NHits[it]; ++i) {
            std::cout << " " << std::setw(4) << HoTArrs[it][i].index;
          }
          std::cout << std::endl;
        }
#endif
      } else {
        add_hit(it, -2, layer_of_hits.layer_id());
        Chi2(it, 0, 0) = 9999.9f;
      }
    }
    dprint("Best candidate = " << BestCandidate(iseed));
  }
}


//==============================================================================
// BestCandidate
//==============================================================================

template<int nseeds, int ncands>
int MkFinderFV<nseeds, ncands>::BestCandidate(int offset) const
{
  int best = -1;
  int nh = 0;
  float chi2 = 9999.9f;
  int it = offset*ncands;
  for (int i = 0; i < ncands; ++i, ++it) {
    if (CandIdx[it] >= 0) {
      const auto nnh = NFoundHits[it];
      const auto nc2 = Chi2[it];
      if (nnh > nh || (nnh == nh && nc2 < chi2)) {
        best = it;
        nh = nnh;
        chi2 = nc2;
      }
    }
  }
  return best;
}
#ifdef INSTANTIATE_FV
template class MkFinderFV<NN/8, 8>;
#endif

} // end namespace mkfit
