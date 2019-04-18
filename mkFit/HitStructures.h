#ifndef HitStructures_H
#define HitStructures_H

#include "../Config.h"
#include "Hit.h"
#include "Track.h"
#include "TrackerInfo.h"
//#define DEBUG
#include "Debug.h"

#include <array>
#include <tbb/tbb.h>

namespace mkfit {

typedef tbb::concurrent_vector<TripletIdx> TripletIdxConVec;

// for each layer
//   Config::nEtaBin vectors of hits, resized to large enough N
//   filled with corresponding hits
//   sorted in phi
//   and with a corresponding BinInfoSomething
//
// My gut feeling is that we could have them all in one place and just store
// indices into small enough eta-phi bins (oh, and eta is a horrible separation variable)
// 10*50000*64/1024/1024
//    30.51757812500000000000
// Or make them smaller ... if we could use short indices that might make a difference.
//
// Need to add eta/phi ... or something significantly better to the Hit struct
//
// At least ... fast eta/phi estimators.
//
// Oh, and we should use radix sort. Could one vectorize this?

// Need a good "array of pods" class with aligned alloc and automatic growth.
// For now just implement the no-resize / no-destroy basics in the BoH.

typedef std::pair<uint16_t, uint16_t> PhiBinInfo_t;

typedef std::vector<PhiBinInfo_t> vecPhiBinInfo_t;

typedef std::vector<vecPhiBinInfo_t> vecvecPhiBinInfo_t;

//==============================================================================

inline bool sortHitsByPhiMT(const Hit& h1, const Hit& h2)
{
  return std::atan2(h1.position()[1],h1.position()[0])<std::atan2(h2.position()[1],h2.position()[0]);
}

inline bool sortTrksByPhiMT(const Track& t1, const Track& t2)
{
  return t1.momPhi() < t2.momPhi();
}

//==============================================================================
//==============================================================================

// Note: the same code is used for barrel and endcap. In barrel the longitudinal
// bins are in Z and in endcap they are in R -- here this coordinate is called Q

class LayerOfHits
{
public:
  const LayerInfo          *m_layer_info = 0;
  Hit                      *m_hits = 0;
  vecvecPhiBinInfo_t        m_phi_bin_infos;
  std::vector<float>        m_hit_phis;
  std::vector<float>        m_hit_qs;

  float m_qmin, m_qmax, m_fq;
  int   m_nq = 0;
  int   m_capacity = 0;

  int   layer_id()  const { return m_layer_info->m_layer_id;    }
  bool  is_barrel() const { return m_layer_info->is_barrel();   }
  bool  is_endcap() const { return ! m_layer_info->is_barrel(); }
  int   bin_index(int q, int p) const { return q*Config::m_nphi + p; }

  PhiBinInfo_t operator[](int i) const {
    int q = i/Config::m_nphi;
    int p = i % Config::m_nphi;
    return m_phi_bin_infos[q][p];
  }

  bool  is_within_z_limits(float z) const { return m_layer_info->is_within_z_limits(z); }
  bool  is_within_r_limits(float r) const { return m_layer_info->is_within_r_limits(r); }

  WSR_Result is_within_z_sensitive_region(float z, float dz) const
  { return m_layer_info->is_within_z_sensitive_region(z, dz); }

  WSR_Result is_within_r_sensitive_region(float r, float dr) const
  { return m_layer_info->is_within_r_sensitive_region(r, dr); }

  float min_dphi() const { return m_layer_info->m_select_min_dphi; }
  float max_dphi() const { return m_layer_info->m_select_max_dphi; }
  float min_dq()   const { return m_layer_info->m_select_min_dq;   }
  float max_dq()   const { return m_layer_info->m_select_max_dq;   }

  // Adding flag for mono/stereo
  bool is_stereo_lyr() const 
  { return  m_layer_info->is_stereo_lyr(); }
 
  // Adding info on sub-detector
  bool is_pixb_lyr() const 
  { return  m_layer_info->is_pixb_lyr(); }
  bool is_pixe_lyr() const 
  { return  m_layer_info->is_pixe_lyr(); }
  bool is_pix_lyr() const 
  { return  m_layer_info->is_pix_lyr(); }
  bool is_tib_lyr() const 
  { return  m_layer_info->is_tib_lyr(); }
  bool is_tob_lyr() const 
  { return  m_layer_info->is_tob_lyr(); }
  bool is_tid_lyr() const 
  { return  m_layer_info->is_tid_lyr(); }
  bool is_tec_lyr() const 
  { return  m_layer_info->is_tec_lyr(); }

  // Adding hit selection limits dynamic factors
  float qf_treg() const { return m_layer_info->m_qf_treg; }
  float phif_treg() const { return m_layer_info->m_phif_treg; }
  float phif_lpt_brl() const { return m_layer_info->m_phif_lpt_brl; }
  float phif_lpt_treg() const { return m_layer_info->m_phif_lpt_treg; }
  float phif_lpt_ec() const { return m_layer_info->m_phif_lpt_ec; }

  // Testing bin filling
  static constexpr float m_fphi     = Config::m_nphi / Config::TwoPI;
  static constexpr int   m_phi_mask = 0x7f;
  static constexpr int   m_phi_bits = 7;
  static constexpr float m_fphi_fine     =  1024 / Config::TwoPI;
  static constexpr int   m_phi_mask_fine = 0x3ff;
  static constexpr int   m_phi_bits_fine = 10;//can't be more than 16
  static constexpr int   m_phi_bits_shift = m_phi_bits_fine - m_phi_bits;
  static constexpr int   m_phi_fine_mask = ~((1 << m_phi_bits_shift) - 1);

protected:

  void setup_bins(float qmin, float qmax, float dq);

  void alloc_hits(int size)
  {
    m_hits = (Hit*) _mm_malloc(sizeof(Hit) * size, 64);
    m_capacity = size;
    for (int ihit = 0; ihit < m_capacity; ihit++){m_hits[ihit] = Hit();} 
    if (Config::usePhiQArrays)
    {
      m_hit_phis.resize(size);
      m_hit_qs  .resize(size);
    }
  }

  void free_hits()
  {
    _mm_free(m_hits);
  }

  void set_phi_bin(int q_bin, int phi_bin, uint16_t &hit_count, uint16_t &hits_in_bin)
  {
    m_phi_bin_infos[q_bin][phi_bin] = { hit_count, hit_count + hits_in_bin };
    hit_count  += hits_in_bin;
    hits_in_bin = 0;
  }

  void empty_phi_bins(int q_bin, int phi_bin_1, int phi_bin_2, uint16_t hit_count)
  {
    for (int pb = phi_bin_1; pb < phi_bin_2; ++pb)
    {
      m_phi_bin_infos[q_bin][pb] = { hit_count, hit_count };
    }
  }

  void empty_q_bins(int q_bin_1, int q_bin_2, uint16_t hit_count)
  {
    for (int qb = q_bin_1; qb < q_bin_2; ++qb)
    {
      empty_phi_bins(qb, 0, Config::m_nphi, hit_count);
    }
  }

public:
  LayerOfHits() {}

  ~LayerOfHits()
  {
    free_hits();
  }

  void  SetupLayer(const LayerInfo &li);

  void  Reset() {}

  float NormalizeQ(float q) const { if (q < m_qmin) return m_qmin; if (q > m_qmax) return m_qmax; return q; }

  int   GetQBin(float q)    const { return (q - m_qmin) * m_fq; }

  int   GetQBinChecked(float q) const { int qb = GetQBin(q); if (qb < 0) qb = 0; else if (qb >= m_nq) qb = m_nq - 1; return qb; }

  // if you don't pass phi in (-pi, +pi), mask away the upper bits using m_phi_mask or use the Checked version.
  int   GetPhiBinFine(float phi) const { return std::floor(m_fphi_fine * (phi + Config::PI)); }
  int   GetPhiBin(float phi) const { return GetPhiBinFine(phi)>>m_phi_bits_shift; }

  int   GetPhiBinChecked(float phi) const { return GetPhiBin(phi) & m_phi_mask; }

  const vecPhiBinInfo_t& GetVecPhiBinInfo(float q) const { return m_phi_bin_infos[GetQBin(q)]; }

  void  SuckInHits(const HitVec &hitv);

  void  SelectHitIndices(float q, float phi, float dq, float dphi, std::vector<int>& idcs, bool isForSeeding=false, bool dump=false);

  void  PrintBins();
};

//==============================================================================

class EventOfHits
{
public:
  std::vector<LayerOfHits> m_layers_of_hits;
  int                      m_n_layers;

public:
  EventOfHits(TrackerInfo &trk_inf);

  void Reset()
  {
    for (auto &i : m_layers_of_hits)
    {
      i.Reset();
    }
  }

  void SuckInHits(int layer, const HitVec &hitv)
  {
    m_layers_of_hits[layer].SuckInHits(hitv);
  }
};



//==============================================================================
// CombinedCandidate and EventOfCombinedCandidates
//==============================================================================

// This inheritance sucks but not doing it will require more changes.

class CombCandidate : public std::vector<Track>
{
public:
  enum SeedState_e { Dormant = 0, Finding, Finished };

  Track        m_best_short_cand;
  SeedState_e  m_state           = Dormant;
  int          m_last_seed_layer = -1;
  unsigned int m_seed_type = 0;

  void MergeCandsAndBestShortOne(bool update_score, bool sort_cands);
};


class EventOfCombCandidates
{
public:
  std::vector<CombCandidate> m_candidates;

  int     m_capacity;
  int     m_size;

public:
  EventOfCombCandidates(int size=0) :
    m_candidates(),
    m_capacity  (0),
    m_size      (0)
  {
    Reset(size);
  }

  void Reset(int new_capacity)
  {
    for (int s = 0; s < m_size; ++s)
    {
      m_candidates[s].clear();
    }

    if (new_capacity > m_capacity)
    {
      m_candidates.resize(new_capacity);

      for (int s = m_capacity; s < new_capacity; ++s)
      {
        m_candidates[s].reserve(Config::maxCandsPerSeed); //we should never exceed this
        m_candidates[s].m_best_short_cand.setCandScore( getScoreWorstPossible() );
      }

      m_capacity = new_capacity;
    }

    m_size = 0;
  }

  CombCandidate& operator[](int i) { return m_candidates[i]; }

  void InsertSeed(const Track& seed)
  {
    assert (m_size < m_capacity);

    m_candidates[m_size].push_back(seed);
    m_candidates[m_size].m_state           = CombCandidate::Dormant;
    m_candidates[m_size].m_last_seed_layer = seed.getLastHitLyr();
    m_candidates[m_size].m_seed_type = seed.getSeedTypeForRanking();
    Track &cand = m_candidates[m_size].back();
    cand.setSeedTypeForRanking(seed.getSeedTypeForRanking());
    cand.setCandScore         (getScoreCand(seed));
    ++m_size;
  }

  void InsertTrack(const Track& track, int seed_index)
  {
    assert (seed_index < m_size);

    m_candidates[seed_index].push_back(track);
  }
};

} // end namespace mkfit
#endif
