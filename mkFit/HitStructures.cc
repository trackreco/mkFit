#include "HitStructures.h"

#include "Event.h"

#include "Ice/IceRevisitedRadix.h"

void LayerOfHits::SetupLayer(float qmin, float qmax, float dq, int layer, bool is_barrel)
{
  // XXXX MT: Could have n_phi bins per layer, too.

  assert (m_nq == 0 && "SetupLayer() already called.");

  m_layer_id  = layer;
  m_is_barrel = is_barrel;

  // Realign qmin/max on "global" dq boundaries.
  // This might not be the best strategy, esp. for large bins.
  // Could have instead:
  // a) extra space on each side;
  // b) pass in number of bins.

  float nmin = std::floor(qmin / dq);
  float nmax = std::ceil (qmax / dq);
  m_qmin = dq * nmin;
  m_qmax = dq * nmax;
  m_fq   = 1.0f / dq; // qbin = (q_hit - m_qmin) * m_fq;
  m_nq   = nmax - nmin;

  m_phi_bin_infos.resize(m_nq);
  for (int i = 0; i < m_nq; ++i) m_phi_bin_infos[i].resize(Config::m_nphi);
}

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

  const int size = hitv.size();

  if (m_capacity < size)
  {
    free_hits();
    alloc_hits(1.02 * size);
  }

#ifndef LOH_USE_PHI_Q_ARRAYS
  std::vector<float>        m_hit_phis(size);
#endif

  struct HitInfo
  {
    float phi;
    float q;
    int   qbin;
  };

  std::vector<HitInfo> ha(size);
  std::vector<int>     qc(m_nq, 0);
  int nqh = m_nq / 2;
  {
    int i = 0;
    for (auto const &h : hitv)
    {
      HitInfo &hi = ha[i];
      hi.phi  = h.phi();
      hi.q    = m_is_barrel ? h.z() : h.r();
      hi.qbin = std::max(std::min(static_cast<int>((hi.q - m_qmin) * m_fq), m_nq - 1), 0);
      m_hit_phis[i] = hi.phi + 6.3f * (hi.qbin - nqh);
      ++qc[hi.qbin];
      ++i;
    }
  }

  RadixSort sort;
  sort.Sort(&m_hit_phis[0], size);

  int curr_q_bin   = 0;
  int curr_phi_bin = 0;
  int hits_in_bin  = 0;
  int hit_count    = 0;

  for (int i = 0; i < size; ++i)
  {
    int j = sort.GetRanks()[i];

    // XXXX MT: Endcap has special check - try to get rid of this!
    // Also, WTF ... this brings in holes as pos i is not filled.
    // If this stays I need i_offset variable.
    if ( ! m_is_barrel && (hitv[j].r() > m_qmax || hitv[j].r() < m_qmin))
    {
      printf("LayerOfHits::SuckInHits WARNING hit out of r boundary of disk\n"
             "  layer %d hit %d hit_r %f limits (%f, %f)\n",
             m_layer_id, j, hitv[j].r(), m_qmin, m_qmax);
      // Figure out of this needs to stay ... and fix it
      // --m_size;
      // continue;
    }

    // Could fix the mis-sorts. Set ha size to size + 1 and fake last entry to avoid ifs.

    memcpy(&m_hits[i], &hitv[j], sizeof(Hit));
#ifdef LOH_USE_PHI_Q_ARRAYS
    m_hit_phis[i] = ha[j].phi;
    m_hit_qs  [i] = ha[j].q;
#endif

    // Fill the bin info

    if (ha[j].qbin != curr_q_bin)
    {
      set_phi_bin(curr_q_bin, curr_phi_bin, hit_count, hits_in_bin);

      empty_phi_bins(curr_q_bin, curr_phi_bin + 1, Config::m_nphi, hit_count);

      empty_q_bins(curr_q_bin + 1, ha[j].qbin, hit_count);

      curr_q_bin = ha[j].qbin;
      curr_phi_bin = 0;
    }

    int phi_bin = GetPhiBin(ha[j].phi);

    if (phi_bin > curr_phi_bin)
    {
      set_phi_bin(curr_q_bin, curr_phi_bin, hit_count, hits_in_bin);

      empty_phi_bins(curr_q_bin, curr_phi_bin + 1, phi_bin, hit_count);

      curr_phi_bin = phi_bin;
    }

    ++hits_in_bin;
  }

  set_phi_bin(curr_q_bin, curr_phi_bin, hit_count, hits_in_bin);

  empty_phi_bins(curr_q_bin, curr_phi_bin + 1, Config::m_nphi, hit_count);

  empty_q_bins(curr_q_bin + 1, m_nq, hit_count);

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
    // XXXX MT: Here I reuse Config::m_max_dz also for endcap ... those will become
    // layer dependant and get stored in TrackerInfo (or copied into LayerOfHits).
    if (std::abs(dq)   > Config::m_max_dz)   dq   = Config::m_max_dz;
    if (std::abs(dphi) > Config::m_max_dphi) dphi = Config::m_max_dphi;
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

      for (int hi = m_phi_bin_infos[qi][pb].first; hi < m_phi_bin_infos[qi][pb].second; ++hi)
      {
        // Here could enforce some furhter selection on hits
#ifdef LOH_USE_PHI_Q_ARRAYS
        float ddq   = std::abs(q   - m_hit_qs[hi]);
        float ddphi = std::abs(phi - m_hit_phis[hi]);
        if (ddphi > Config::PI) ddphi = Config::TwoPI - ddphi;

        if (dump)
          printf("     SHI %3d %4d %4d %5d  %6.3f %6.3f %6.4f %7.5f   %s\n",
                 qi, pi, pb, hi,
                 m_hit_qs[hi], m_hit_phis[hi], ddq, ddphi,
                 (ddq < dq && ddphi < dphi) ? "PASS" : "FAIL");

        if (ddq < dq && ddphi < dphi)
#endif
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
    printf("%c bin %d\n", m_is_barrel ? 'Z' : 'R', qb);
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

EventOfHits::EventOfHits(int n_layers) :
  m_layers_of_hits(n_layers),
  m_n_layers(n_layers)
{
  for (int i = 0; i < n_layers; ++i)
  {
    assert("Butchering out old setup via Config ... for cyl-cow and cms." && 0);
    //if (Config::endcapTest) m_layers_of_hits[i].SetupDisk(Config::cmsDiskMinRs[i], Config::cmsDiskMaxRs[i], Config::g_disk_dr[i]);
    //else m_layers_of_hits[i].SetupLayer(-Config::g_layer_zwidth[i], Config::g_layer_zwidth[i], Config::g_layer_dz[i]);
  }
}

EventOfHits::EventOfHits(TrackerInfo &trk_inf) :
  m_layers_of_hits(trk_inf.m_layers.size()),
  m_n_layers(trk_inf.m_layers.size())
{
  for (auto &li : trk_inf.m_layers)
  {
    // XXXX MT bin width has to come from somewhere ... had arrays in Config.
    // This is also different for strip detectors + the hack
    // with mono / stereo regions in CMS endcap layers.
    // I'm just hardocding this to 1 cm for now ...

    float bin_width = 1.0f;

    if (li.is_barrel())
    {
      // printf("EventOfHits::EventOfHits setting up layer %2d as barrel\n", li.m_layer_id);

      m_layers_of_hits[li.m_layer_id].SetupLayer(li.m_zmin, li.m_zmax, bin_width,
                                                 li.m_layer_id, true);
    }
    else
    {
      // printf("EventOfHits::EventOfHits setting up layer %2d as endcap\n", li.m_layer_id);

      m_layers_of_hits[li.m_layer_id].SetupLayer(li.m_rin, li.m_rout, bin_width,
                                                 li.m_layer_id, false);
    }
  }
}
