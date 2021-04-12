#include "HitStructures.h"

#include "Event.h"

#include "SteeringParams.h" // just needs iteration-params

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
}

void LayerOfHits::SetupLayer(const LayerInfo &li)
{
  // Note, LayerInfo::m_q_bin ==>  > 0 - bin width, < 0 - number of bins

  assert (m_layer_info == nullptr && "SetupLayer() already called.");

  m_layer_info = &li;

  m_is_barrel = m_layer_info->is_barrel();

  if (m_is_barrel) setup_bins(li.m_zmin, li.m_zmax, li.m_q_bin);
  else             setup_bins(li.m_rin,  li.m_rout, li.m_q_bin);
}

//==============================================================================

/*
void detect_q_min_max(const HitVec &hitv)
{
  float dq = 0.5;
  // should know this from geom.
  //m_qmin =  1000;
  //m_qmax = -1000;
  for (auto const &h : hitv)
  {
    if (h.q() < m_qmin) m_qmin = h.q();
    if (h.q() > m_qmax) m_qmax = h.q();
  }
  printf("LoH::SuckInHits qmin=%f, qmax=%f", m_qmin, m_qmax);
  float nmin = std::floor(m_qmin / dq);
  float nmax = std::ceil (m_qmax / dq);
  m_qmin = dq * nmin;
  m_qmax = dq * nmax;
  int nq  = nmax - nmin;
  int nqh = nq / 2;
  m_fq = 1.0f / dq; // qbin = (qhit - m_qmin) * m_fq;
  printf(" -> qmin=%f, qmax=%f, nq=%d, fq=%f\n", m_qmin, m_qmax, nq, m_fq);
}
*/

void LayerOfHits::SuckInHits(const HitVec &hitv)
{
  assert (m_nq > 0 && "SetupLayer() was not called.");

  const int size = hitv.size();

  m_ext_hits  = & hitv;

#ifdef COPY_SORTED_HITS
  if (m_capacity < size)
  {
    free_hits();
    alloc_hits(1.02 * size);
  }
#endif

  if (Config::usePhiQArrays)
  {
    m_hit_phis.resize(size);
    m_hit_qs.resize(size);
    m_hit_infos.resize(size);
  }
  m_qphifines.resize(size);

  for (int i = 0; i < size; ++i)
  {
    const Hit &h = hitv[i];

    HitInfo hi = { h.phi(), m_is_barrel ? h.z() : h.r() };

    m_qphifines[i] = GetPhiBinFine(hi.phi) + (GetQBinChecked(hi.q) << 16);

    if (Config::usePhiQArrays)
    {
      m_hit_infos[i] = hi;
    }
  }

  operator delete [] (m_hit_ranks);
  {
    RadixSort sort;
    sort.Sort(&m_qphifines[0], size, RADIX_UNSIGNED);
    m_hit_ranks = sort.RelinquishRanks();
  }

  int curr_qphi = -1;
  empty_q_bins(0, m_nq, 0);

  for (int i = 0; i < size; ++i)
  {
    int j = m_hit_ranks[i];

#ifdef COPY_SORTED_HITS
    memcpy(&m_hits[i], &hitv[j], sizeof(Hit));
#endif

    if (Config::usePhiQArrays)
    {
      m_hit_phis[i] = m_hit_infos[j].phi;
      m_hit_qs  [i] = m_hit_infos[j].q;
    }

    // Combined q-phi bin with fine part masked off
    const int jqphi   = m_qphifines[j] & m_phi_fine_xmask;

    const int phi_bin = (jqphi & m_phi_mask_fine) >> m_phi_bits_shift;
    const int q_bin   = jqphi >> 16;

    // Fill the bin info
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

//==============================================================================


void LayerOfHits::BeginRegistrationOfHits(const HitVec &hitv)
{
  assert (m_nq > 0 && "SetupLayer() was not called.");

  m_ext_hits = &hitv;

  m_hit_infos.clear();
  m_qphifines.clear();
  m_ext_idcs.clear();
  m_min_ext_idx = std::numeric_limits<int>::max();
  m_max_ext_idx = std::numeric_limits<int>::min();
}

void LayerOfHits::RegisterHit(int idx)
{
  const Hit &h = (*m_ext_hits)[idx];

  m_ext_idcs.push_back(idx);
  m_min_ext_idx = std::min(m_min_ext_idx, idx);
  m_max_ext_idx = std::max(m_max_ext_idx, idx);

  HitInfo hi = { h.phi(), m_is_barrel ? h.z() : h.r() };

  m_qphifines.push_back( GetPhiBinFine(hi.phi) + (GetQBinChecked(hi.q) << 16) );

  if (Config::usePhiQArrays)
  {
    m_hit_infos.emplace_back(hi);
  }
}

void LayerOfHits::EndRegistrationOfHits(bool build_original_to_internal_map)
{
  const int size = m_ext_idcs.size();
  if (size == 0) return;

  // radix
  operator delete [] (m_hit_ranks);
  {
    RadixSort sort;
    sort.Sort(&m_qphifines[0], size, RADIX_UNSIGNED);
    m_hit_ranks = sort.RelinquishRanks();
  }

  // copy q/phi

#ifdef COPY_SORTED_HITS
  if (m_capacity < size)
  {
    free_hits();
    alloc_hits(1.02 * size);
  }
#endif

  if (Config::usePhiQArrays)
  {
    m_hit_phis.resize(size);
    m_hit_qs.resize(size);
  }

  int curr_qphi = -1;
  empty_q_bins(0, m_nq, 0);

  for (int i = 0; i < size; ++i)
  {
    int j = m_hit_ranks[i];  // index in intermediate
    int k = m_ext_idcs[j];   // index in external hit_vec

#ifdef COPY_SORTED_HITS
    memcpy(&m_hits[i], &hitv[k], sizeof(Hit));
#endif

    if (Config::usePhiQArrays)
    {
      m_hit_phis[i] = m_hit_infos[j].phi;
      m_hit_qs  [i] = m_hit_infos[j].q;
    }

    // Combined q-phi bin with fine part masked off
    const int jqphi = m_qphifines[j] & m_phi_fine_xmask;

    const int phi_bin = (jqphi & m_phi_mask_fine) >> m_phi_bits_shift;
    const int q_bin   = jqphi >> 16;

    // Fill the bin info
    if (jqphi != curr_qphi)
    {
      m_phi_bin_infos[q_bin][phi_bin] = {i, i};
      curr_qphi = jqphi;
    }

    m_phi_bin_infos[q_bin][phi_bin].second++;

    // m_hit_ranks[i] will never be used again ... use it to point to external/original index,
    // as it does in standalone case.
    m_hit_ranks[i] = k;
  }

  if (build_original_to_internal_map)
  {
    if (m_max_ext_idx - m_min_ext_idx + 1 > 8*size)
    {
      // If this happens we might:
      // a) Use external indices for everything. -- *** We are now. ***
      // b) Build these maps for seeding layers only.
      // c) Have a flag in hit-on-track that tells us if the hit index has been remapped,
      //    essentially, if it is a seed hit. This might not be too stupid anyway.
      //    Hmmh, this could be index < -256, or something like that.

      printf("LayerOfHits::EndRegistrationOfHits() original_to_internal index map vector is largish: size=%d, map_vector_size=%d\n",
             size, m_max_ext_idx - m_min_ext_idx + 1);
    }

    m_ext_idcs.resize(m_max_ext_idx - m_min_ext_idx + 1);
    for (int i = 0; i < size; ++i)
    {
      m_ext_idcs[m_hit_ranks[i] - m_min_ext_idx] = i;
    }
  }

  // For Matti: We can release m_hit_infos and m_qphifines -- and realloc on next BeginInput.
}

//==============================================================================


/*
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
*/

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
// TrackCand
//==============================================================================

Track TrackCand::exportTrack() const
{
  // printf("TrackCand::exportTrack label=%5d, total_hits=%2d, overlaps=%2d\n", label(),
  //        nTotalHits(), nOverlapHits_);

  Track res(*this);
  res.resizeHits(nTotalHits(), nFoundHits());
  res.setNOverlapHits(nOverlapHits());

  int nh = nTotalHits();
  int ch = lastHitIdx_;
  while (--nh >= 0)
  {
    HoTNode& hot_node = m_comb_candidate->m_hots[ch];

    res.setHitIdxAtPos(nh, hot_node.m_hot);

    // printf("  nh=%2d, ch=%d, idx=%d lyr=%d prev_idx=%d\n",
    //        nh, ch, hot_node.m_hot.index, hot_node.m_hot.layer, hot_node.m_prev_idx);

    ch       = hot_node.m_prev_idx;
  }

  return res;
}


//==============================================================================
// EventOfHits
//==============================================================================

EventOfHits::EventOfHits(const TrackerInfo &trk_inf) :
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

void CombCandidate::ImportSeed(const Track& seed)
{
  emplace_back(TrackCand(seed, this));

  m_state           = CombCandidate::Dormant;
  m_last_seed_layer = seed.getLastHitLyr();
  m_seed_type       = seed.getSeedTypeForRanking();
  m_seed_algo       = seed.algoint();
  m_seed_label      = seed.label();

  TrackCand &cand = back();

  // printf("Importing pt=%f eta=%f, lastCcIndex=%d\n", cand.pT(), cand.momEta(), cand.lastCcIndex());

  for (const HitOnTrack* hp = seed.BeginHitsOnTrack(); hp != seed.EndHitsOnTrack(); ++hp)
  {
    // printf(" hit idx=%d lyr=%d\n", hp->index, hp->layer);
    cand.addHitIdx(hp->index, hp->layer, 0.0f);
  }

  cand.setSeedTypeForRanking(m_seed_type);
  cand.setScore             (getScoreCand(cand));
}

void CombCandidate::MergeCandsAndBestShortOne(const IterationParams& params, bool update_score, bool sort_cands)
{
  TrackCand *best_short = m_best_short_cand.combCandidate() ? & m_best_short_cand : nullptr;

  if ( ! empty())
  {
    if (update_score)
    {
      for (auto &c : *this) c.setScore( getScoreCand(c) );
      if (best_short) best_short->setScore( getScoreCand(*best_short) );
    }
    if (sort_cands)
    {
      std::sort(begin(), end(), sortByScoreTrackCand);
    }

    if (best_short && best_short->score() > back().score())
    {
      auto ci = begin();
      while (ci->score() > best_short->score()) ++ci;

      if ((int) size() >= params.maxCandsPerSeed) pop_back();

      // To print out what has been replaced -- remove when done with short track handling.
      // if (ci == finalcands.begin()) {
      //   printf("FindTracksStd -- Replacing best cand (%f) with short one (%f) in final sorting\n",
      //          front().score(), best_short->score());
      // }

      insert(ci, *best_short);
    }

  }
  else if (best_short)
  {
    push_back( *best_short );
  }

  if (best_short) best_short->resetShortTrack();

  // assert(capacity() == (size_t)Config::maxCandsPerSeed);
}

} // end namespace mkfit
