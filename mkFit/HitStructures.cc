#include "HitStructures.h"

#include "Event.h"

#include "Ice/IceRevisitedRadix.h"

namespace mkfit {

void LayerOfHits::setup_bins(float qmin, float qmax, float dq)
{
  // Define layer with min/max and number of bins along q.

  // XXXX MT: Could have n_phi bins per layer, too.

  if (dq < 0)
  {
    m_nq = (int) -dq;
    m_qmin = qmin;
    m_qmax = qmax;
  }
  else
  {
    float extent = qmax - qmin;
    m_nq = std::ceil(extent / dq);
    float extra  = 0.5f * (m_nq * dq - extent);
    m_qmin = qmin - extra;
    m_qmax = qmax + extra;
  }
  m_fq = m_nq / (qmax - qmin); // qbin = (q_hit - m_qmin) * m_fq;

  m_phi_bin_infos.resize(m_nq);
  for (int i = 0; i < m_nq; ++i) m_phi_bin_infos[i].resize(Config::m_nphi);
}

void LayerOfHits::SetupLayer(const LayerInfo &li)
{
  // Note, LayerInfo::m_q_bin ==>  > 0 - bin width, < 0 - number of bins

  assert (m_layer_info == nullptr && "SetupLayer() already called.");

  m_layer_info = &li;

  if (is_barrel()) setup_bins(li.m_zmin, li.m_zmax, li.m_q_bin);
  else             setup_bins(li.m_rin,  li.m_rout, li.m_q_bin);
}

//==============================================================================

void LayerOfHits::SuckInHits(const HitVec &hitv)
{
  // This is now in SetupLayer()
  // // should be layer dependant
  // float dq = 0.5;
  // // should know this from geom.
  // //m_qmin =  1000;
  // //m_qmax = -1000;
  // for (auto const &h : hitv)
  // {
  //   if (h.q() < m_qmin) m_qmin = h.q();
  //   if (h.q() > m_qmax) m_qmax = h.q();
  // }
  // printf("LoH::SuckInHits qmin=%f, qmax=%f", m_qmin, m_qmax);
  // float nmin = std::floor(m_qmin / dq);
  // float nmax = std::ceil (m_qmax / dq);
  // m_qmin = dq * nmin;
  // m_qmax = dq * nmax;
  // int nq  = nmax - nmin;
  // int nqh = nq / 2;
  // m_fq = 1.0f / dq; // qbin = (qhit - m_qmin) * m_fq;
  // printf(" -> qmin=%f, qmax=%f, nq=%d, fq=%f\n", m_qmin, m_qmax, nq, m_fq);

  assert (m_nq > 0 && "SetupLayer() was not called.");

  const int  size   = hitv.size();
  const bool is_brl = is_barrel();

  if (m_capacity < size)
  {
    free_hits();
    alloc_hits(1.02 * size);
  }

  if (!Config::usePhiQArrays)
  {
    m_hit_phis.resize(size);
  }

  struct HitInfo
  {
    float phi;
    float q;
    uint16_t   qbin;
    uint16_t   phibin;
  };

  std::vector<HitInfo> ha(size);
  std::vector<udword>     hit_qphiFines(size);
  
  {
    for (int i =0; i < size; ++i)
    {
      auto const& h = hitv[i];
      
      HitInfo &hi = ha[i];
      // N.1.a For phi in [-pi, pi): squashPhiMinimal(h.phi()); Apparently atan2 can round the wrong way.
      hi.phi  = h.phi();
      hi.q    = is_brl ? h.z() : h.r();
      hi.qbin = std::max(std::min(static_cast<int>((hi.q - m_qmin) * m_fq), m_nq - 1), 0);
      // N.1.b ha[j].phi is returned by atan2 and can be rounded the wrong way, resulting in bin -1 or Config::m_nphi
      hi.phibin = GetPhiBin(hi.phi) & m_phi_mask;
      hit_qphiFines[i] = GetPhiBinFine(hi.phi) + (hi.qbin << 16);
    }//i < size
  }

  RadixSort sort;
  sort.Sort(&hit_qphiFines[0], size);

  int curr_qphi    = -1;
  empty_q_bins(0, m_nq, 0);

  for (uint16_t i = 0; i < size; ++i)
  {
    int j = sort.GetRanks()[i];

    // XXXX MT: Endcap has special check - try to get rid of this!
    // Also, WTF ... this brings in holes as pos i is not filled.
    // If this stays I need i_offset variable.
    /* UNCOMMENT FOR DEBUGS??
    if ( ! is_brl && (hitv[j].r() > m_qmax || hitv[j].r() < m_qmin))
    {
      printf("LayerOfHits::SuckInHits WARNING hit out of r boundary of disk\n"
             "  layer %d hit %d hit_r %f limits (%f, %f)\n",
             layer_id(), j, hitv[j].r(), m_qmin, m_qmax);
      // Figure out of this needs to stay ... and fix it
      // --m_size;
      // continue;
    }
    */

    // Could fix the mis-sorts. Set ha size to size + 1 and fake last entry to avoid ifs.

    memcpy(&m_hits[i], &hitv[j], sizeof(Hit));
    if (Config::usePhiQArrays)
    {
      m_hit_phis[i] = ha[j].phi;
      m_hit_qs  [i] = ha[j].q;
    }

    const int phi_bin = ha[j].phibin;
    const int q_bin = ha[j].qbin;

    // Fill the bin info
    const int jqphi = hit_qphiFines[j] & m_phi_fine_mask;
    if (jqphi != curr_qphi)
    {
      m_phi_bin_infos[q_bin][phi_bin] = {i, i};
      curr_qphi = jqphi;
    }

    m_phi_bin_infos[q_bin][phi_bin].second++;
  }


  // Check for mis-sorts due to lost precision (not really important).
  // float phi_prev = 0;
  // int   bin_prev = -1;
  // int   prev_err_idx = -1;
  // for (int i = 0; i < size; ++i)
  // {
  //   int j = sort.GetRanks()[i];
  //   float phi  = ha[j].phi;
  //   int   qbin = ha[j].qbin;
  //   if (qbin == bin_prev && phi < phi_prev)
  //   {
  //     //printf("  Offset error: %5d %5d %10f %10f %5d %8f\n", i, j, phi, phi_prev, ha[j].qbin, hitv[j].q());
  //     if (prev_err_idx == i - 1)
  //       printf("DOUBLE Offset error: %5d %5d %10f %10f %5d %8f\n", i, j, phi, phi_prev, ha[j].qbin, hitv[j].q());
  //     prev_err_idx = i;
  //   }
  //   phi_prev = phi;
  //   bin_prev = qbin;
  // }

  // Print first couple
  // for(int i = 0; i < 20; ++i)
  // {
  //   int j = sort.GetRanks()[i];
  //
  //   printf("%3d %3d %8f %5d %8f\n", i, j, ha[j].phi, ha[j].qbin, hitv[j].q());
  // }
}

void LayerOfHits::SelectHitIndices(float q, float phi, float dq, float dphi, std::vector<int>& idcs, bool isForSeeding, bool dump)
{
  // Sanitizes q, dq and dphi. phi is expected to be in -pi, pi.

  // Make sure how phi bins work beyond -pi, +pi.
  // for (float p = -8; p <= 8; p += 0.05)
  // {
  //   int pb = GetPhiBin(p);
  //   printf("%5.2f %4d %4d\n", p, pb, pb & m_phi_mask);
  // }

  if ( ! isForSeeding) // seeding has set cuts for dq and dphi
  {
    // XXXX MT: min search windows not enforced here.
    dq   = std::min(std::abs(dq),   max_dq());
    dphi = std::min(std::abs(dphi), max_dphi());
  }
  
  int qb1 = GetQBinChecked(q - dq);
  int qb2 = GetQBinChecked(q + dq) + 1;
  int pb1 = GetPhiBin(phi - dphi);
  int pb2 = GetPhiBin(phi + dphi) + 1;

  // int extra = 2;
  // qb1 -= 2; if (qb < 0) qb = 0;
  // qb2 += 2; if (qb >= m_nq) qb = m_nq;

  if (dump)
    printf("LayerOfHits::SelectHitIndices %6.3f %6.3f %6.4f %7.5f %3d %3d %4d %4d\n",
           q, phi, dq, dphi, qb1, qb2, pb1, pb2);

  // This should be input argument, well ... it will be Matriplex op, or sth. // KPM -- it is now! used for seeding
  for (int qi = qb1; qi < qb2; ++qi)
  {
    for (int pi = pb1; pi < pb2; ++pi)
    {
      int pb = pi & m_phi_mask;

      for (uint16_t hi = m_phi_bin_infos[qi][pb].first; hi < m_phi_bin_infos[qi][pb].second; ++hi)
      {
        // Here could enforce some furhter selection on hits
	if (Config::usePhiQArrays)
	{
	  float ddq   = std::abs(q   - m_hit_qs[hi]);
	  float ddphi = std::abs(phi - m_hit_phis[hi]);
	  if (ddphi > Config::PI) ddphi = Config::TwoPI - ddphi;
	  
	  if (dump)
	    printf("     SHI %3d %4d %4d %5d  %6.3f %6.3f %6.4f %7.5f   %s\n",
		   qi, pi, pb, hi,
		   m_hit_qs[hi], m_hit_phis[hi], ddq, ddphi,
		   (ddq < dq && ddphi < dphi) ? "PASS" : "FAIL");
	  
	  if (ddq < dq && ddphi < dphi)
	  {
	    idcs.push_back(hi);
	  }
	}
	else // do not use phi-q arrays
	{
	  idcs.push_back(hi);
	}
      }
    }
  }
}

void LayerOfHits::PrintBins()
{
  for (int qb = 0; qb < m_nq; ++qb)
  {
    printf("%c bin %d\n", is_barrel() ? 'Z' : 'R', qb);
    for (int pb = 0; pb < Config::m_nphi; ++pb)
    {
      if (pb % 8 == 0)
        printf(" Phi %4d: ", pb);
      printf("%5d,%4d   %s",
             m_phi_bin_infos[qb][pb].first, m_phi_bin_infos[qb][pb].second,
             ((pb + 1) % 8 == 0) ? "\n" : "");
    }
  }
}


//==============================================================================
// EventOfHits
//==============================================================================

EventOfHits::EventOfHits(TrackerInfo &trk_inf) :
  m_layers_of_hits(trk_inf.m_layers.size()),
  m_n_layers(trk_inf.m_layers.size())
{
  for (auto &li : trk_inf.m_layers)
  {
    m_layers_of_hits[li.m_layer_id].SetupLayer(li);
  }
}


//==============================================================================
// CombCandidate
//==============================================================================

void CombCandidate::MergeCandsAndBestShortOne(bool update_score, bool sort_cands)
{
  std::vector<Track> &finalcands = *this;
  Track              &best_short = m_best_short_cand;

  if ( ! finalcands.empty())
  {
    if (update_score)
    {
      for (auto &c : finalcands) c.setCandScore( getScoreCand(c) );
    }
    if (sort_cands)
    {
      std::sort(finalcands.begin(), finalcands.end(), sortByScoreCand);
    }

    if (best_short.getCandScore() > finalcands.back().getCandScore())
    {
      auto ci = finalcands.begin();
      while (ci->getCandScore() > best_short.getCandScore()) ++ci;

      if (finalcands.size() > Config::maxCandsPerSeed)  finalcands.pop_back();

      // To print out what has been replaced -- remove when done with short track handling.
      /*
        if (ci == finalcands.begin())
        {
        printf("FindTracksStd -- Replacing best cand (%d) with short one (%d) in final sorting for seed index=%d\n",
                     finalcands.front().getCandScore(), best_short.getCandScore(), iseed);
        }
      */

      finalcands.insert(ci, best_short);
    }

  }
  else if (best_short.getCandScore() > getScoreWorstPossible())
  {
    finalcands.push_back( best_short );
  }

  best_short.setCandScore( getScoreWorstPossible() );
}

} // end namespace mkfit
